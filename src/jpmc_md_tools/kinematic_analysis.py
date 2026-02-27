import numpy as np
from typing import Optional
import logging

logger = logging.getLogger(__name__)

def COM_trj(positions_oxygens: np.array, positions_h1s: np.array, positions_h2s: np.array, velocities_oxygens: np.array, velocities_h1s: np.array, velocities_h2s: np.array, data: dict, hydrogen_type: int, oxygen_type: int) -> tuple:
    '''
    Calculate the center of mass (COM) trajectory for a system of water molecules.
    Parameters:
        positions_oxygens: A numpy array of shape (M, N_water_molecules, 3) representing the positions of oxygen atoms.
        positions_h1s: A numpy array of shape (M, N_water_molecules, 3) representing the positions of the first hydrogen atoms.
        positions_h2s: A numpy array of shape (M, N_water_molecules, 3) representing the positions of the second hydrogen atoms.
        velocities_oxygens: A numpy array of shape (M, N_water_molecules, 3) representing the velocities of oxygen atoms.
        velocities_h1s: A numpy array of shape (M, N_water_molecules, 3) representing the velocities of the first hydrogen atoms.
        velocities_h2s: A numpy array of shape (M, N_water_molecules, 3) representing the velocities of the second hydrogen atoms.
        data: A dictionary containing the atomic data, including masses and atom types.
        hydrogen_type: An integer representing the atom type for hydrogen in the data dictionary.
        oxygen_type: An integer representing the atom type for oxygen in the data dictionary.
    Returns:
        A tuple containing:
            - positions_COM: A numpy array of shape (M, N_water_molecules, 3) representing the positions of the center of mass for each water molecule.
            - velocities_COM: A numpy array of shape (M, N_water_molecules, 3) representing the velocities of the center of mass for each water molecule.
    '''
    
    hydrogen_mass=data["Masses"][hydrogen_type]
    oxygen_mass=data["Masses"][oxygen_type]

    positions_COM=(positions_oxygens*oxygen_mass+positions_h1s*hydrogen_mass+positions_h2s*hydrogen_mass)/(oxygen_mass+2*hydrogen_mass)
    velocities_COM=(velocities_oxygens*oxygen_mass+velocities_h1s*hydrogen_mass+velocities_h2s*hydrogen_mass)/(oxygen_mass+2*hydrogen_mass)

    return positions_COM, velocities_COM

def arrange_trj_data_by_molecules(data: dict, trj: np.array, oxygen_type: int, hydrogen_type: int, global2local: Optional[dict] = None) -> tuple:
    '''
    Arrange trajectory data by molecules for a system of water molecules.
    Parameters:
        data: A dictionary containing the atomic data, including atom types and molecule IDs.
        trj: A numpy array of shape (M, N_atoms, 4) representing the trajectory data, where M is the number of time steps, N_atoms is the total number of atoms, and the last dimension contains the time, x, y, z coordinates.
        oxygen_type: An integer representing the atom type for oxygen in the data dictionary.
        hydrogen_type: An integer representing the atom type for hydrogen in the data dictionary.
        global2local: An optional dictionary mapping global atom IDs to local atom IDs, used to filter the trajectory data to only include certain atoms.
    Returns:
        A tuple containing:
            - oxygens: A numpy array of shape (M, N_water_molecules, 3) representing the positions of oxygen atoms for each water molecule.
            - h1s: A numpy array of shape (M, N_water_molecules, 3) representing the positions of the first hydrogen atoms for each water molecule.
            - h2s: A numpy array of shape (M, N_water_molecules, 3) representing the positions of the second hydrogen atoms for
    '''
    
    # Create mappings
    atom_2_molecule = {}
    OXYGENS = {}
    H1S = {}
    H2S = {}
    N_water_molecules = 0
    local_molecule_id=0
    molecule_id_map={}

    for atom in data["Atoms"]:
        atom_2_molecule[atom] = data["Atoms"][atom]["molecule ID"]
        if data["Atoms"][atom]["atom type"] == oxygen_type:
            N_water_molecules += 1
            OXYGENS[atom_2_molecule[atom]] = atom
        if data["Atoms"][atom]["atom type"] == hydrogen_type:
            if atom_2_molecule[atom] not in H1S:
                H1S[atom_2_molecule[atom]] = atom
            elif atom != H1S[atom_2_molecule[atom]]:
                H2S[atom_2_molecule[atom]] = atom

    # remake OXYGENS, H1S, and H2S to only include the molecules that have at least one atom in the global2local directory
    if not global2local==None:
        new_molecule_ID=0
        old2new_moleculeID_map={}
        newOXYGENS={}
        newH1S={}
        newH2S={}
        for molecule in OXYGENS:
            oxygen = OXYGENS[molecule]
            h1 = H1S[molecule]
            h2 = H2S[molecule]
            if (oxygen in global2local) and (h1 in global2local) and (h2 in global2local):
                new_molecule_ID+=1
                old2new_moleculeID_map[molecule]=new_molecule_ID
                newOXYGENS[new_molecule_ID]=global2local[OXYGENS[molecule]]
                newH1S[new_molecule_ID]=global2local[H1S[molecule]]
                newH2S[new_molecule_ID]=global2local[H2S[molecule]]
        N_water_molecules=new_molecule_ID
        OXYGENS=newOXYGENS
        H1S=newH1S
        H2S=newH2S

    # Get dimensions
    M = trj.shape[0]

    # Initialize arrays
    oxygens = np.zeros((M, N_water_molecules, 3))
    h1s = np.zeros((M, N_water_molecules, 3))
    h2s = np.zeros((M, N_water_molecules, 3))

    # Populate arrays based on molecule indexing
    for mol_id in range(1, N_water_molecules + 1):  # Assuming molecule IDs start at 1
        if mol_id in OXYGENS:
            oxygens[:, mol_id - 1, :] = np.double(trj[:, OXYGENS[mol_id] - 1, 1:])
        if mol_id in H1S:
            h1s[:, mol_id - 1, :] = np.double(trj[:, H1S[mol_id] - 1, 1:])
        if mol_id in H2S:
            h2s[:, mol_id - 1, :] = np.double(trj[:, H2S[mol_id] - 1, 1:])
        
    return oxygens, h1s, h2s

def compute_local_basis_unit_vectors(data: dict, trj: np.array, oxygen_type: int, hydrogen_type: int, box_size: float, mode: Optional[str]="debug", global2local: Optional[dict]=None) -> tuple:
    '''
    Compute local basis unit vectors for a system of water molecules based on their positions in the trajectory data.
    Parameters:
        data: A dictionary containing the atomic data, including atom types and molecule IDs.
        trj: A numpy array of shape (M, N_atoms, 4) representing the trajectory data, where M is the number of time steps, N_atoms is the total number of atoms, and the last dimension contains the time, x, y, z coordinates.
        oxygen_type: An integer representing the atom type for oxygen in the data dictionary.
        hydrogen_type: An integer representing the atom type for hydrogen in the data dictionary.
        box_size: A float representing the size of the simulation box, used for applying periodic boundary conditions when calculating relative positions.
        mode: An optional string parameter that can be set to "debug" to enable additional print statements and plots for debugging purposes.
        global2local: An optional dictionary mapping global atom IDs to local atom IDs, used to filter the trajectory data to only include certain atoms.
    Returns:
        A tuple containing:
            - a: A numpy array of shape (M, N_water_molecules, 3) representing the first local basis unit vector for each water molecule at each time step.
            - b: A numpy array of shape (M, N_water_molecules, 3) representing the second local basis unit vector for each water molecule at each time step.
            - c: A numpy array of shape (M, N_water_molecules, 3) representing the third local basis unit vector for each water molecule at each time step.
            - oxygens: A numpy array of shape (M, N_water_molecules, 3) representing the positions of oxygen atoms for each water molecule
    '''


    oxygens, h1s, h2s = arrange_trj_data_by_molecules(data,trj,oxygen_type,hydrogen_type,global2local=global2local)

    r_h1_rel = h1s - oxygens - np.round((h1s - oxygens) / box_size) * box_size
    r_h2_rel = h2s - oxygens - np.round((h2s - oxygens) / box_size) * box_size

    logger.debug("maximum oh1 length: %f", np.max(np.linalg.norm(r_h1_rel,axis=2)))
    logger.debug("maximum oh2 length: %f", np.max(np.linalg.norm(r_h2_rel,axis=2)))

    # FIND NEW BASIS VECTORS

    a = r_h1_rel+r_h2_rel
    a = np.divide(a, np.linalg.norm(a,axis=2,keepdims=True), where=np.linalg.norm(a,axis=2,keepdims=True) != 0)  # Avoid division by zero (set to 0 where norm is zero)

    b = np.cross(r_h1_rel,a,axis=2)
    b = np.divide(b, np.linalg.norm(b,axis=2,keepdims=True), where=np.linalg.norm(b,axis=2,keepdims=True) != 0)  # Avoid division by zero (set to 0 where norm is zero)

    c = np.cross(a,b,axis=2)
    c = np.divide(c, np.linalg.norm(c,axis=2,keepdims=True), where=np.linalg.norm(c,axis=2,keepdims=True) != 0)  # Avoid division by zero (set to 0 where norm is zero)

    # CHECKS
    
    test_molecule=0

    logger.debug("a:")
    logger.debug(a[0,test_molecule,:])
    logger.debug("norm(a):")
    logger.debug(np.linalg.norm(a[0,test_molecule,:]))
    logger.debug("b:")
    logger.debug(b[0,test_molecule,:])
    logger.debug("norm(b):")
    logger.debug(np.linalg.norm(b[0,test_molecule,:]))
    logger.debug("c:")
    logger.debug(c[0,test_molecule,:])
    logger.debug("norm(c):")
    logger.debug(np.linalg.norm(c[0,test_molecule,:]))
    logger.debug("a dot b:")
    logger.debug(np.dot(a[0,test_molecule,:],b[0,test_molecule,:]))
    logger.debug("b dot c:")
    logger.debug(np.dot(b[0,test_molecule,:],c[0,test_molecule,:]))
    logger.debug("c dot a:")
    logger.debug(np.dot(a[0,test_molecule,:],c[0,test_molecule,:]))

    return a, b, c, oxygens

def compute_de_dt(a: np.array, b: np.array, c: np.array, dt: float, style: Optional[str]="central difference") -> tuple:
    '''
    Compute the time derivatives of the local basis unit vectors for a system of water molecules.
    Parameters:
        a: A numpy array of shape (M, N_water_molecules, 3) representing the first local basis unit vector for each water molecule at each time step.
        b: A numpy array of shape (M, N_water_molecules, 3) representing the second local basis unit vector for each water molecule at each time step.
        c: A numpy array of shape (M, N_water_molecules, 3) representing the third local basis unit vector for each water molecule at each time step.
        dt: A float representing the time step between consecutive frames in the trajectory data.
        style: An optional string parameter that can be set to "backward difference", "forward difference", or "central difference" to specify the finite difference method used for computing the derivatives. The default is "central difference".
    Returns:
        A tuple containing:
            - da_dt: A numpy array of shape (M, N_water_molecules, 3) representing the time derivative of the first local basis unit vector for each water molecule at each time step.
            - db_dt: A numpy array of shape (M, N_water_molecules, 3) representing the time derivative of the second local basis unit vector for each water molecule at each time step.
            - dc_dt: A numpy array of shape (M, N_water_molecules, 3) representing the time derivative of the third local basis unit vector for each water molecule at each time step.
    '''
    
    M=a.shape[0]
    da_dt=np.zeros(a.shape)
    db_dt=np.zeros(b.shape)
    dc_dt=np.zeros(c.shape)

    if style=="backward difference":
        da_dt[1:,:,:]=(a[1:,:,:]-a[:M-1,:,:])/dt
        db_dt[1:,:,:]=(b[1:,:,:]-b[:M-1,:,:])/dt
        dc_dt[1:,:,:]=(c[1:,:,:]-c[:M-1,:,:])/dt

    elif style=="forward difference":
        da_dt[:M-1,:,:]=(a[1:,:,:]-a[:M-1,:,:])/dt
        db_dt[:M-1,:,:]=(b[1:,:,:]-b[:M-1,:,:])/dt
        dc_dt[:M-1,:,:]=(c[1:,:,:]-c[:M-1,:,:])/dt

    elif style=="central difference":
        da_dt[1:M-1,:,:]=(a[2:,:,:]-a[:M-2,:,:])/(2*dt)
        db_dt[1:M-1,:,:]=(b[2:,:,:]-b[:M-2,:,:])/(2*dt)
        dc_dt[1:M-1,:,:]=(c[2:,:,:]-c[:M-2,:,:])/(2*dt)

    return da_dt, db_dt, dc_dt

def compute_angular_velocity_from_basis_vectors(a: np.array, b: np.array, c: np.array, da_dt: np.array, db_dt: np.array, dc_dt: np.array, style: Optional[str]="central difference") -> np.array:
    '''
    Compute the angular velocity of the local basis vectors for a system of water molecules.
    Parameters:
        a: A numpy array of shape (M, N_water_molecules, 3) representing the first local basis unit vector for each water molecule at each time step.
        b: A numpy array of shape (M, N_water_molecules, 3) representing the second local basis unit vector for each water molecule at each time step.
        c: A numpy array of shape (M, N_water_molecules, 3) representing the third local basis unit vector for each water molecule at each time step.
        da_dt: A numpy array of shape (M, N_water_molecules, 3) representing the time derivative of the first local basis unit vector for each water molecule at each time step.
        db_dt: A numpy array of shape (M, N_water_molecules, 3) representing the time derivative of the second local basis unit vector for each water molecule at each time step.
        dc_dt: A numpy array of shape (M, N_water_molecules, 3) representing the time derivative of the third local basis unit vector for each water molecule at each time step.
        style: An optional string parameter that can be set to "backward difference", "forward difference", or "central difference" to specify the finite difference method used for computing the derivatives. The default is "central difference". This parameter is used to determine how to slice the output array to match the time steps for which the derivatives were computed.
    Returns:
        A numpy array of shape (M, N_water_molecules, 3) representing the angular velocity of the local basis vectors for each water molecule at each time step. The last dimension contains the three components of the angular velocity vector. The time steps included in the output array depend on the finite difference method used for computing the derivatives, as specified by the style parameter. For "backward difference", the first time step is excluded; for "forward difference", the last time step is excluded; for "central difference", the first and last time steps are excluded.
    '''
    
    M=a.shape[0]

    angular_velocities=np.zeros(a.shape)

    angular_velocities[:,:,0]=np.sum(db_dt*c,axis=2)
    angular_velocities[:,:,1]=np.sum(dc_dt*a,axis=2)
    angular_velocities[:,:,2]=np.sum(da_dt*b,axis=2)

    if style=="backward difference":
        angular_velocities=angular_velocities[1:,:,:]

    elif style=="forward difference":
        angular_velocities=angular_velocities[:M-1,:,:]

    elif style=="central difference":
        angular_velocities=angular_velocities[1:M-1,:,:]

    return angular_velocities

def unwrap_trj(oxygens: np.array, h1s: np.array, h2s: np.array, data: dict) -> tuple:
    '''
    Apply periodic boundary conditions to unwrap the trajectory data for a system of water molecules.
    Parameters:
        oxygens: A numpy array of shape (M, N_water_molecules, 3) representing the positions of oxygen atoms for each water molecule at each time step.
        h1s: A numpy array of shape (M, N_water_molecules, 3) representing the positions of the first hydrogen atoms for each water molecule at each time step.
        h2s: A numpy array of shape (M, N_water_molecules, 3) representing the positions of the second hydrogen atoms for each water molecule at each time step.
        data: A dictionary containing the atomic data, including the dimensions of the simulation box (xlo, xhi, ylo, yhi, zlo, zhi) used for applying periodic boundary conditions.
    Returns:
        A tuple containing:
            - oxygens: A numpy array of shape (M, N_water_molecules, 3) representing the unwrapped positions of oxygen atoms for each water molecule at each time step. The positions are adjusted to account for periodic boundary conditions, ensuring that the trajectory data reflects the true movement of the molecules across the simulation box boundaries.
            - h1s: A numpy array of shape (M, N_water_molecules, 3) representing the unwrapped positions of the first hydrogen atoms for each water molecule at each time step, adjusted for periodic boundary conditions.
            - h2s: A numpy array of shape (M, N_water_molecules, 3) representing the unwrapped positions of the second hydrogen atoms for each water molecule at each time step, adjusted for periodic boundary conditions.
    '''

    box_size=[np.double(data["xhi"])-np.double(data["xlo"]),np.double(data["yhi"])-np.double(data["ylo"]),np.double(data["zhi"])-np.double(data["zlo"])]

    h1s = h1s - np.round((h1s - oxygens) / box_size) * box_size
    h2s = h2s - np.round((h2s - oxygens) / box_size) * box_size

    return oxygens,h1s,h2s

def compute_internal_velocities(local_atom_v: np.array, local_COM_v: np.array, r_local: np.array, omega: np.array) -> np.array:
    """
    Compute the internal velocities of atoms in a molecule by removing the contributions from the center of mass motion and the rotational motion.
    Parameters:
        local_atom_v: A numpy array of shape (T, M, 3) representing the velocities of individual atoms in the local basis for each molecule at each time step.
        local_COM_v: A numpy array of shape (T, M, 3) representing the velocities of the center of mass for each molecule at each time step.
        r_local: A numpy array of shape (T, M, 3) representing the positions of individual atoms relative to the center of mass in the local basis for each molecule at each time step.
        omega: A numpy array of shape (T, M, 3) representing the angular velocity of the local basis vectors for each molecule at each time step.
    Returns:
        A numpy array of shape (T, M, 3) representing the internal velocities of the atoms in the local basis for each molecule at each time step, with the contributions from the center of mass motion and the rotational motion removed. The internal velocity is calculated as:
        internal_velocity = local_atom_v - local_COM_v - (omega × r_local)
        where "×" denotes the cross product. This calculation effectively isolates the vibrational motion of the atoms within the molecule by removing the translational motion of the center of mass and the rotational motion of the molecule as a whole.
    """
    # Cross product of angular velocity with local position: ω × r
    # omega: shape (T, M, 3)
    # r_local: shape (T, M, 3)

    rotational_contrib = np.cross(omega, r_local)  # shape: (T, M, 3)

    return local_atom_v - local_COM_v - rotational_contrib

def project_onto_local_ref_frame(v: np.array, a: np.array, b: np.array, c: np.array) -> np.array:
    '''
    Project the velocities of atoms onto the local reference frame defined by the basis vectors a, b, and c for each molecule at each time step.
    Parameters:
        v: A numpy array of shape (T, M, 3) representing the velocities of individual atoms in the global reference frame for each molecule at each time step.
        a: A numpy array of shape (T, M, 3) representing the first local basis unit vector for each molecule at each time step.
        b: A numpy array of shape (T, M, 3) representing the second local basis unit vector for each molecule at each time step.
        c: A numpy array of shape (T, M, 3) representing the third local basis unit vector for each molecule at each time step.
    Returns:
        A numpy array of shape (T, M, 3) representing the velocities of individual atoms projected onto the local reference frame defined by the basis vectors a, b, and c for each molecule at each time step. The projection is calculated by taking the dot product of the velocity vector with each of the basis vectors, effectively transforming the velocity from the global reference frame to the local reference frame defined by the orientation of the molecule. The resulting array contains the components of the velocity in the directions of the local basis vectors a, b, and c, which can be used for further analysis of the internal dynamics of the molecule.
    '''
    
    
    if v.shape[2]>3:
        print("WARNING: v has shape "+str(v.shape))
        print("Using v[:,:,"+str(v.shape[2]-3)+":]")
        v=v[:,:,v.shape[2]-3]

    # Each of A, B, C has shape (M, N, 3)
    R = np.stack([a, b, c], axis=-2)  # Shape: (M, N, 3, 3)

    # v: shape (M, N, 3)
    # R: shape (M, N, 3, 3)
    local_v = np.einsum('mnij,mnj->mni', R, v)

    return local_v





