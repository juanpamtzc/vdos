from jpmc_md_tools.io_tools import readTRJFile, readDatFile
from jpmc_md_tools.kinematic_analysis import (
    compute_local_basis_unit_vectors,
    arrange_trj_data_by_molecules,
    unwrap_trj,
    COM_trj,
    project_onto_local_ref_frame,
    compute_de_dt,
    compute_angular_velocity_from_basis_vectors,
    compute_internal_velocities,
    compute_vdos_fft,
)
from jpmc_md_tools.spectral_analysis import compute_vdos_from_velocity_fft
import numpy as np
import argparse
import os

if __name__ == "__main__":
    # --- Set up Command Line Arguments ---
    parser = argparse.ArgumentParser(description="Compute normalized VDOS spectra from LAMMPS trajectories.")
    
    # Required Arguments
    parser.add_argument("--data_dir", type=str, required=True, help="Path to the directory containing trajectory files.")
    parser.add_argument("--filename_prefix", type=str, required=True, help="Prefix for the files (e.g., '0' for '0Cpsqm_...').")
    
    # Optional Arguments with sensible defaults
    parser.add_argument("--charge_sign", type=str, choices=["bulk", "positive", "negative"], default="bulk", help="Which interfacial region to analyze.")
    parser.add_argument("--dt", type=float, default=0.5, help="Timestep between frames in femtoseconds.")
    parser.add_argument("--oxygen_type", type=int, default=2, help="LAMMPS atom type ID for Oxygen.")
    parser.add_argument("--hydrogen_type", type=int, default=3, help="LAMMPS atom type ID for Hydrogen.")
    parser.add_argument("--interfacial_limit", type=float, default=5.0, help="Distance limit for the interface mask.")
    
    args = parser.parse_args()

    # --- Use the Arguments ---
    # Ensure the path ends with a slash
    path = os.path.join(args.data_dir, "") 
    filename = args.filename_prefix
    charge_sign = args.charge_sign
    dt = args.dt
    oxygen_type = args.oxygen_type
    hydrogen_type = args.hydrogen_type
    interfacial_limit = args.interfacial_limit

    # Construct file names dynamically
    position_fn = f"{path}{filename}Cpsqm_1.0M_KCl_{charge_sign}_pos_VDOS.lammpstrj"
    velocity_fn = f"{path}{filename}Cpsqm_1.0M_KCl_{charge_sign}_vel_VDOS.lammpstrj"
    data_fn = f"{path}{filename}Cpsqm_1.0M_KCl_eq_1.dat"

    print(f"--- Starting VDOS Analysis ---")
    print(f"Data Dir: {path}")
    print(f"Region: {charge_sign}")

    # Load data
    data = readDatFile(data_fn)
    positions, local2global, global2local = readTRJFile(position_fn, True, None)
    velocities, _, _ = readTRJFile(velocity_fn, True, None)

    box_size = [
        np.double(data["xhi"]) - np.double(data["xlo"]),
        np.double(data["yhi"]) - np.double(data["ylo"]),
        np.double(data["zhi"]) - np.double(data["zlo"]),
    ]

    n_timesteps = positions.shape[0]
    n_atoms = data["# atoms"]
    



    # Compute local orthonormal basis vectors (per molecule, per timestep)
    a, b, c, _ = compute_local_basis_unit_vectors(
        data, positions, oxygen_type, hydrogen_type, box_size, mode="debug", global2local=global2local
    )

    # Split positions/velocities by atom types in molecules
    positions_oxygens, positions_h1s, positions_h2s = arrange_trj_data_by_molecules(
        data, positions, oxygen_type, hydrogen_type, global2local=global2local
    )

    velocities_oxygens, velocities_h1s, velocities_h2s = arrange_trj_data_by_molecules(
        data, velocities, oxygen_type, hydrogen_type, global2local=global2local
    )

    # Remove PBC effects (unwrap positions)
    positions_oxygens, positions_h1s, positions_h2s = unwrap_trj(
        positions_oxygens, positions_h1s, positions_h2s, data
    )

    # Compute center-of-mass positions and velocities
    positions_COM, velocities_COM = COM_trj(
        positions_oxygens, positions_h1s, positions_h2s,
        velocities_oxygens, velocities_h1s, velocities_h2s,
        data, hydrogen_type, oxygen_type
    )

    # interfacial masks
    if charge_sign == "negative":
        interfacial_mask = positions_COM[0,:,2] <= interfacial_limit
    elif charge_sign == "positive":
        interfacial_mask = positions_COM[0,:,2] >= 100 - interfacial_limit
    elif charge_sign == "bulk":
        interfacial_mask = (positions_COM[0,:,2] >= 45.0) & (positions_COM[0,:,2] <= 55.0)
            
    # Compute positions relative to COM (in global frame)
    rel_positions_oxygens = positions_oxygens - positions_COM
    rel_positions_h1s = positions_h1s - positions_COM
    rel_positions_h2s = positions_h2s - positions_COM

    # Project relative positions into local frame
    local_positions_oxygens = project_onto_local_ref_frame(rel_positions_oxygens, a, b, c)
    local_positions_h1s = project_onto_local_ref_frame(rel_positions_h1s, a, b, c)
    local_positions_h2s = project_onto_local_ref_frame(rel_positions_h2s, a, b, c)

    # Project velocities into local frame
    local_velocities_oxygens = project_onto_local_ref_frame(velocities_oxygens, a, b, c)
    local_velocities_h1s = project_onto_local_ref_frame(velocities_h1s, a, b, c)
    local_velocities_h2s = project_onto_local_ref_frame(velocities_h2s, a, b, c)
    local_velocities_COM = project_onto_local_ref_frame(velocities_COM, a, b, c)

    # Compute angular velocities from rotating basis
    da_dt, db_dt, dc_dt = compute_de_dt(a, b, c, dt, style="central difference")
    angular_velocities = compute_angular_velocity_from_basis_vectors(
        a, b, c, da_dt, db_dt, dc_dt, style="central difference"
    )  # shape: (T, M, 3)

    # VDOS Computations
    suffix = f".npy"

    #compute_vdos_from_velocity_fft(velocities_COM[:,interfacial_mask,:], dt, output_vdos_file=path + filename + '_translational_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(local_velocities_COM[:,interfacial_mask,0], dt, output_vdos_file=path + filename + '_translational_a_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(local_velocities_COM[:,interfacial_mask,1], dt, output_vdos_file=path + filename + '_translational_b_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(local_velocities_COM[:,interfacial_mask,2], dt, output_vdos_file=path + filename + '_translational_c_vdos_from_spectrum_' + charge_sign + suffix)

    #compute_vdos_from_velocity_fft(velocities_oxygens[:,interfacial_mask,:], dt, output_vdos_file=path + filename + '_oxygen_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(local_velocities_oxygens[:,interfacial_mask,0], dt, output_vdos_file=path + filename + '_oxygen_a_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(local_velocities_oxygens[:,interfacial_mask,1], dt, output_vdos_file=path + filename + '_oxygen_b_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(local_velocities_oxygens[:,interfacial_mask,2], dt, output_vdos_file=path + filename + '_oxygen_c_vdos_from_spectrum_' + charge_sign + suffix)
    #compute_vdos_from_velocity_fft(velocities_h1s[:,interfacial_mask,:], dt, output_vdos_file=path + filename + '_h1s_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(local_velocities_h1s[:,interfacial_mask,0], dt, output_vdos_file=path + filename + '_h1s_a_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(local_velocities_h1s[:,interfacial_mask,1], dt, output_vdos_file=path + filename + '_h1s_b_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(local_velocities_h1s[:,interfacial_mask,2], dt, output_vdos_file=path + filename + '_h1s_c_vdos_from_spectrum_' + charge_sign + suffix)
    #compute_vdos_from_velocity_fft(velocities_h2s[:,interfacial_mask,:], dt, output_vdos_file=path + filename + '_h2s_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(local_velocities_h2s[:,interfacial_mask,0], dt, output_vdos_file=path + filename + '_h2s_a_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(local_velocities_h2s[:,interfacial_mask,1], dt, output_vdos_file=path + filename + '_h2s_b_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(local_velocities_h2s[:,interfacial_mask,2], dt, output_vdos_file=path + filename + '_h2s_c_vdos_from_spectrum_' + charge_sign + suffix)

    freqs, vdos_data = compute_vdos_fft(angular_velocities[:,interfacial_mask,:], dt, output_vdos_file=path + filename + '_angular_vdos_from_spectrum_' + charge_sign + suffix)

    internal_velocities_oxygens = compute_internal_velocities(
        local_velocities_oxygens[1:-1], local_velocities_COM[1:-1], local_positions_oxygens[1:-1], angular_velocities
    )
    internal_velocities_h1s = compute_internal_velocities(
        local_velocities_h1s[1:-1], local_velocities_COM[1:-1], local_positions_h1s[1:-1], angular_velocities
    )
    internal_velocities_h2s = compute_internal_velocities(
        local_velocities_h2s[1:-1], local_velocities_COM[1:-1], local_positions_h2s[1:-1], angular_velocities
    )

    internal_velocities_all = np.concatenate([
        np.double(data["Masses"][oxygen_type]) * internal_velocities_oxygens[:,interfacial_mask,:],
        np.double(data["Masses"][hydrogen_type]) * internal_velocities_h1s[:,interfacial_mask,:],
        np.double(data["Masses"][hydrogen_type]) * internal_velocities_h2s[:,interfacial_mask,:],
    ], axis=1)

    compute_vdos_from_velocity_fft(internal_velocities_all, dt, output_vdos_file=path + filename + '_internal_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(internal_velocities_all[:,:,0], dt, output_vdos_file=path + filename + '_internal_a_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(internal_velocities_all[:,:,1], dt, output_vdos_file=path + filename + '_internal_b_vdos_from_spectrum_' + charge_sign + suffix)
    compute_vdos_from_velocity_fft(internal_velocities_all[:,:,2], dt, output_vdos_file=path + filename + '_internal_c_vdos_from_spectrum_' + charge_sign + suffix)
        
    print("--- VDOS Analysis Complete ---")
