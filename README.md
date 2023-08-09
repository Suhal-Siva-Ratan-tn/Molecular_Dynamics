# Molecular_Dynamics(MD)

## Introduction
The primary goal of this project is to simulate the interactions and movements of particles (atoms, molecules) over time using the Molecular Dynamics approach. The simulation calculates forces between particles, integrates their equations of motion, and produces trajectories that capture the system's behavior. This repository provides the code in FORTRAN, C, and CUDA C languages. The provided FORTRAN and C files for Molecular Dynamics implementation also offer compatibility with OPENACC mode.

## Usage

Run the following commands based on the file

To compile the `md.c` file, use the following command to include the math library: 

```bash
gcc md.c -lm
