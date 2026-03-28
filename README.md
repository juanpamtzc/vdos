# High-Performance LAMMPS VDOS Analysis Pipeline

[![ci](https://github.com/juanpamtzc/vdos/actions/workflows/ci.yml/badge.svg)](https://github.com/juanpamtzc/vdos/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11-blue.svg)
![LAMMPS](https://img.shields.io/badge/LAMMPS-MD-red.svg)
![HPC](https://img.shields.io/badge/HPC-SLURM-orange)

## Overview
This repository contains a modular, high-performance Python pipeline for extracting, processing, and analyzing large-scale Molecular Dynamics (MD) trajectories generated via LAMMPS. 

Specifically, this codebase computes the Vibrational Density of States (VDOS) of water by evaluating the mass-weighted velocity power spectral density (PSD):

$$g(\tilde{\nu}) = \frac{2}{k_B T} \sum_{i=1}^{N} \sum_{\alpha=1}^{d} m_i \lim_{\tau \to \infty} \frac{1}{2\tau} \mathbb{E} \left[ \left| \int_{-\tau}^{\tau} v_i^\alpha(t) e^{-i 2\pi \nu t} dt \right|^2 \right]$$

The VDOS gives spectral information of the fundamental dynamical modes of water (i.e., the frequencies at which different movements occur).

To give a more detailed glimpse into the nature of the various dynamical modes, this pipeline allows the VDOS to be decomposed per degree of freedom, differentiating between translational, rotational, and internal degrees of freedom (as well as their directional components).

## Algorithmic Optimization: The Wiener-Khinchin Theorem
A naive computation of the VDOS requires computing the Velocity Autocorrelation Function (VACF), which scales at $O(N^2)$ for trajectory length $N$. 

By applying the **Wiener-Khinchin theorem**, the pipeline computes the VDOS directly from the magnitude squared of the Fast Fourier Transform (FFT) of the velocity trajectories. This reduces the time complexity to **$O(N \log N)$**, resulting in massive performance gains without any loss of mathematical fidelity.

**For the rigorous mathematical proof, see:** [`docs/VDOS_Mathematical_Proof.pdf`](docs/VDOS_Mathematical_Proof.pdf)

## Repository Architecture

```text
├── .github/workflows/       # CI/CD pipelines
├── docs/                    # Mathematical documentation and proofs
├── notebooks/               # Jupyter Notebooks for publication-ready visualizations
├── scripts/                 # Execution scripts for local and cluster environments
│   └── compute_vdos.py      # Main entry point for VDOS calculation
├── src/jpmc_md_tools/       # Core computational logic (The Engine)
│   ├── __init__.py
│   ├── lammps_io.py         # Memory-efficient I/O for LAMMPS formats
│   ├── kinematics.py        # Spatial math, internal velocities, and COM translations
│   └── spectra.py           # FFT logic, signal processing, and VDOS
├── tests/                   # Unit test suite (pytest/unittest)
├── pyproject.toml           # Modern Python packaging configuration
└── requirements.txt         # Dependency management