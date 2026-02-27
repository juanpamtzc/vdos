# High-Performance LAMMPS VDOS Analysis Pipeline

![Python](https://img.shields.io/badge/python-3.9+-blue.svg)
![LAMMPS](https://img.shields.io/badge/LAMMPS-MD-red.svg)
![HPC](https://img.shields.io/badge/HPC-SLURM-orange)

## Overview
This repository contains a modular, high-performance Python pipeline for extracting, processing, and analyzing large-scale Molecular Dynamics (MD) trajectories generated via LAMMPS. 

Specifically, this codebase processes physical systems at $1 \text{ atm}$ (e.g., interfacial water and ion clusters), isolating internal molecular kinematics from center-of-mass (COM) and rigid-body rotational motion. It computes the Vibrational Density of States (VDOS) by computing the mass-weighed velocity power spectrum:

$$\text{VDOS}(\omega) \propto \int_{-\infty}^{\infty} \langle \vec{v}(t) \cdot \vec{v}(0) \rangle e^{-i\omega t} dt$$

This project was built to process multi-gigabyte trajectory datasets efficiently on cluster supercomputers using batch processing, while outputting lightweight data for visualization.

## 📁 Repository Architecture

```text
├── lammps_inputs/           # LAMMPS initial configuration and execution scripts
├── data/
│   ├── raw/                 # Lightweight starting configurations (.dat)
│   └── trajectories/        # [GIT IGNORED] Heavy .lammpstrj binary outputs
├── notebooks/               # Jupyter Notebooks for final publication figure generation
├── scripts/                 # Execution scripts for local and cluster environments
│   ├── 01_parse_trajectories.py
│   ├── 02_compute_vdos_bulk.py
│   ├── 03_compute_vdos_ions.py
│   └── run_lammps_cluster.sh
├── src/                     # Core computational logic (The Engine)
│   ├── lammps_io.py         # Memory-efficient I/O for LAMMPS formats
│   ├── kinematics.py        # Spatial math, internal velocities, and COM translations
│   └── spectra.py           # FFT logic and VDOS computations
└── README.md