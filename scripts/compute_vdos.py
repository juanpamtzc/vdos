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

if __name__ == "__main__":

    trajectory_path="/media/jpmc/easystore/work/EDL_OSC/channel/1M_KCl/trajectories/VDOS_5/"
    data_path="/media/jpmc/easystore/work/EDL_OSC/channel/1M_KCl/trajectories/VDOS_5/"
    output_path="/media/jpmc/easystore/work/EDL_OSC/channel/1M_KCl/trajectories/VDOS_5/"
    
    filenames=["0"]
    charge_signs=["bulk"]
    path="/media/jpmc/easystore/work/EDL_OSC/channel/1M_KCl/trajectories/VDOS_5/"
    oxygen_type=2
    hydrogen_type=3
    dt=0.5
    interfacial_limit=5.0

    for filename in filenames:
        for charge_sign in charge_signs:    
            position_fn=path+filename+"Cpsqm_1.0M_KCl_"+charge_sign+"_pos_VDOS.lammpstrj"
            velocity_fn=path+filename+"Cpsqm_1.0M_KCl_"+charge_sign+"_vel_VDOS.lammpstrj"
            data_fn=path+filename+"Cpsqm_1.0M_KCl_eq_1.dat"

            # Load data
            data = readDatFile(data_fn)
            positions, local2global, global2local = readTRJFile(position_fn, True, None)
            velocities, _, _ = readTRJFile(velocity_fn, True, None)
            #positions = np.load(position_fn)
            #velocities = np.load(velocity_fn)

            box_size = [
                np.double(data["xhi"]) - np.double(data["xlo"]),
                np.double(data["yhi"]) - np.double(data["ylo"]),
                np.double(data["zhi"]) - np.double(data["zlo"]),
            ]

            n_timesteps = positions.shape[0]
            n_atoms=data["# atoms"]
            mid = n_timesteps // 2

            for i in range(2):
                # Select time slice
                if i == 0:
                    pos_chunk = positions[:mid,:,:]
                    vel_chunk = velocities[:mid,:,:]
                else:
                    pos_chunk = positions[mid:,:,:]
                    vel_chunk = velocities[mid:,:,:]

                # Compute local orthonormal basis vectors (per molecule, per timestep)
                a, b, c, _ = compute_local_basis_unit_vectors(
                    data, pos_chunk, oxygen_type, hydrogen_type, box_size, mode="debug", global2local=global2local
                )

                # Split positions/velocities by atom types in molecules
                positions_oxygens, positions_h1s, positions_h2s = arrange_trj_data_by_molecules(
                    data, pos_chunk, oxygen_type, hydrogen_type, global2local=global2local
                )
                velocities_oxygens, velocities_h1s, velocities_h2s = arrange_trj_data_by_molecules(
                    data, vel_chunk, oxygen_type, hydrogen_type, global2local=global2local
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
                if charge_sign=="negative":
                    interfacial_mask = positions_COM[0,:,2]<=interfacial_limit
                elif charge_sign=="positive":
                    interfacial_mask = positions_COM[0,:,2]>=100-interfacial_limit
                elif charge_sign=="bulk":
                    interfacial_mask = (positions_COM[0,:,2]>=45.0) & (positions_COM[0,:,2]<=55.0)
                    
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
                suffix = f"_chunk_{i}.npy"

                compute_vdos_from_velocity_fft(velocities_COM[:,interfacial_mask,:], dt, output_vdos_file=path + filename + '_translational_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(local_velocities_COM[:,interfacial_mask,0], dt, output_vdos_file=path + filename + '_translational_a_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(local_velocities_COM[:,interfacial_mask,1], dt, output_vdos_file=path + filename + '_translational_b_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(local_velocities_COM[:,interfacial_mask,2], dt, output_vdos_file=path + filename + '_translational_c_vdos_from_spectrum' + charge_sign + suffix)

                compute_vdos_from_velocity_fft(velocities_oxygens[:,interfacial_mask,:], dt, output_vdos_file=path + filename + '_oxygen_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(local_velocities_oxygens[:,interfacial_mask,0], dt, output_vdos_file=path + filename + '_oxygen_a_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(local_velocities_oxygens[:,interfacial_mask,1], dt, output_vdos_file=path + filename + '_oxygen_b_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(local_velocities_oxygens[:,interfacial_mask,2], dt, output_vdos_file=path + filename + '_oxygen_c_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(velocities_h1s[:,interfacial_mask,:], dt, output_vdos_file=path + filename + '_h1s_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(local_velocities_h1s[:,interfacial_mask,0], dt, output_vdos_file=path + filename + '_h1s_a_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(local_velocities_h1s[:,interfacial_mask,1], dt, output_vdos_file=path + filename + '_h1s_b_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(local_velocities_h1s[:,interfacial_mask,2], dt, output_vdos_file=path + filename + '_h1s_c_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(velocities_h2s[:,interfacial_mask,:], dt, output_vdos_file=path + filename + '_h2s_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(local_velocities_h2s[:,interfacial_mask,0], dt, output_vdos_file=path + filename + '_h2s_a_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(local_velocities_h2s[:,interfacial_mask,1], dt, output_vdos_file=path + filename + '_h2s_b_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(local_velocities_h2s[:,interfacial_mask,2], dt, output_vdos_file=path + filename + '_h2s_c_vdos_from_spectrum' + charge_sign + suffix)

                freqs, vdos_data = compute_vdos_fft(angular_velocities[:,interfacial_mask,:], dt, output_vdos_file=path + filename + '_angular_vdos_from_spectrum' + charge_sign + suffix)

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

                compute_vdos_from_velocity_fft(internal_velocities_all, dt, output_vdos_file=path + filename + '_internal_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(internal_velocities_all[:,:,0], dt, output_vdos_file=path + filename + '_internal_a_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(internal_velocities_all[:,:,1], dt, output_vdos_file=path + filename + '_internal_b_vdos_from_spectrum' + charge_sign + suffix)
                compute_vdos_from_velocity_fft(internal_velocities_all[:,:,2], dt, output_vdos_file=path + filename + '_internal_c_vdos_from_spectrum' + charge_sign + suffix)