# Molecular_Dynamics(MD)

## Introduction
The primary goal of this project is to simulate the interactions and movements of particles (atoms, molecules) over time using the Molecular Dynamics approach. The simulation calculates forces between particles, integrates their equations of motion, and produces trajectories that capture the system's behavior. This repository provides the code in FORTRAN, C, and CUDA C languages. The provided FORTRAN and C files for Molecular Dynamics implementation also offer compatibility with OPENACC mode.

## Usage

To compile the `md.c` file, use the following command to include the math library: 

```bash
gcc md.c -lm
```

To compile the `md.c` file with OpenACC support, use the following command:

```bash
nvc -acc -Minfo=accel md.c
```

To compile the `md2.f90` file, use the following command:

```bash
gfortran md2.f90
```

To compile the `md2.f90` file with OpenACC support, use the following command:

```bash
nvfortran -acc -Minfo=accel md2.f90
```

To compile the `md.cu` file, use the following command to include cuBLAS library:

```bash
nvcc md.cu -lcublas
```
