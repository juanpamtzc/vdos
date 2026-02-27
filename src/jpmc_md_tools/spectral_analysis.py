import numpy as np
from scipy.fft import fft, fftfreq
from typing import Optional

def compute_vdos_from_velocity_fft(velocities: np.array, dt: np.double, output_vdos_file: Optional[str]='vdos_from_vel_fft.npy') -> tuple:
    """
    Assumes that the given velocities are for only one particle species (e.g., oxygen atoms, water molecules, etc.)

    Computes the vibrational density of states (VDOS) by taking the Fourier transform
    of the velocities directly and calculating the power spectrum.

    Parameters:
        velocities: numpy array of shape [N_t, N_atoms, 3] containing vx, vy, vz for selected atoms.
        dt: time step between frames.
        output_vdos_file: Filename to save the VDOS data (default: 'vdos_from_vel_fft.npy').
    Returns:
        A tuple containing:
            - positive_freqs: numpy array of positive frequencies.
            - positive_power_spectrum: numpy array of corresponding VDOS values (power spectrum).
    Saves:
        VDOS (from velocity FFT) as a numpy array file with the frequency and corresponding VDOS values.
    """

    N_t=velocities.shape[0]
    N_atoms=velocities.shape[1]
    total_velocities_flat = velocities.reshape(N_t, -1)  # Flatten velocity components for all atoms

    # Compute the FFT of the velocities
    velocity_fft = fft(total_velocities_flat, axis=0)

    # Compute the power spectrum (squared magnitude of the FFT)
    power_spectrum = np.abs(velocity_fft)**2

    # Calculate the frequencies
    freqs = fftfreq(N_t, d=dt)

    # Consider only the positive frequencies
    positive_freq_mask = freqs >= 0
    positive_freqs = freqs[positive_freq_mask]
    positive_power_spectrum = power_spectrum[positive_freq_mask].mean(axis=1) # Average over velocity components

    vdos_data = np.column_stack((positive_freqs, positive_power_spectrum))
    np.save(output_vdos_file, vdos_data)

    return positive_freqs, positive_power_spectrum


