# Impurity-Preserved Partition Based Density Matrix Embedding Theory for Local Electronic Excitations

Computational data for the manuscript **"Impurity-Preserved Partition Based Density Matrix Embedding Theory for Local Electronic Excitations"**.

Manuscript tracking number: ct-2026-012516

## Overview

This repository contains the input files, checkpoint files and representative output files for the quantum-chemical calculations reported in the manuscript. The calculations include all-electron reference calculations and DMET calculations with different partitioning schemes for the lanthanide single-ion magnets and Ce(III) luminescent complexes. These calculations are performed using our in-house Python package [LAMP_dmet](https://github.com/ccme-tmc/LAMP_dmet).
The package is built on [PySCF](https://github.com/pyscf/pyscf).

The `embed_sim/` directory contains the development version of LAMP_dmet when these calculations were performed.

The `iao_helper.py` script is an in-house Python extension of the intrinsic and projected atomic orbitals (IAO+PAO) generation module in the [QC-DMET](https://github.com/sebwouters/qc-dmet) package, developed to support lanthanide molecules and open-shell systems.

For the `1Dy/` benchmark set, the all-electron results in both `CAS/` and `PT2/` were taken from our
previous work, DOI:
[10.1021/acs.jctc.5c01336](https://doi.org/10.1021/acs.jctc.5c01336). 
The input files and representative output files for the corresponding calcualtions are also included.



## Directory Layout

- `1Dy/`: Calculations for
  $[\mathrm{Dy}(\mathrm{Cp}^{\mathrm{ttt}})_2]^+$.

    The all-electron results in both `1Dy/PT2/AE` and `1Dy/CAS/AE` were taken from our
    previous work, DOI:[10.1021/acs.jctc.5c01336](https://doi.org/10.1021/acs.jctc.5c01336). 
    The input files and representative output files for the corresponding calculations are also included. The corresponding output files are provided as `1Dy/PT2/AE/jobid*.out`. All PT2 correction energies were summarized in `1Dy/PT2/AE/nevpt2.txt`, and the corrected energies were summarized in `1Dy/PT2/AE/DyCP2_opt.txt`. The subsequent results of two-step SISO calculations were summarized in `1Dy/PT2/AE/DyCP2_mag.txt`.

- `2Dy/`: Calculations for
  $[\mathrm{Dy}\{\mathrm{N}[\mathrm{Si}({}^{i}\mathrm{Pr})_3][\mathrm{Si}({}^{i}\mathrm{Pr})_2\mathrm{C}(\mathrm{CH}_3)=\mathrm{CHCH}_3]\}]^+$.
  
  In `2Dy/PT2/AE`, the all-electron NEVPT2-SOC calculations of 2Dy were performed at the
  same theoretical level as for the other systems. The only procedural difference
  was introduced to reduce the memory demand of the challenging NEVPT2 calculation
  with 854 basis functions.
  First, the PT2 correction energies were computed for states 1–16; the
  corresponding output is provided as the representative output file
  `2Dy/PT2/AE/1_16pt2.out`. The PT2 correction energies for the remaining states
  17–21 were then computed separately. All PT2 correction energies were summarized
  in `2Dy/PT2/AE/Dy_nevpt2.txt`, which was subsequently used as the input for the
  two-step SISO calculation to obtain the final results.
  Most calculations in this repository used the development version of the package
  included in `embed_sim/`. For the all-electron NEVPT2-SOC calculation of 2Dy, however, a newer
  version of [LAMP_dmet](https://github.com/ccme-tmc/LAMP_dmet) was used to enable
  state-specific NEVPT2 calculations. The corresponding package files  are included in `2Dy/PT2/AE/embed_sim_new/`.

- `3Ce/`: Calculations for
  $\mathrm{Ce\text{-}Bp}^{\mathrm{Me}}$.

- `4Ce/`: Calculations for
  $\mathrm{Ce\text{-}O}_2\mathrm{ip}^{\mathrm{tBu}}$.

- `embed_sim/`: In-house Python package for the embedding calculations (the
  version used when these calculations were performed).

For each molecular system, the calculation folders are organized according to the
high-level solvers used. Common subdirectories include:

- `CAS/`: CASSCF-SOC calculations.

- `PT2/`: NEVPT2-SOC calculations.

Within each high-level-solver directory, the calculation folders are further
organized by the use of all-electron or DMET calculations, the impurity defined
in the DMET calculations, and the partitioning scheme used in DMET. Common
subdirectories include:

- `AE/`: All-electron calculations with the corresponding high-level solver.

- `AO_small/`: Atomic orbital based DMET (AO-DMET) calculations with the central
  lanthanide ion as the impurity.

- `AO_big/`: AO-DMET calculations with the lanthanide
  ion and its nearest coordinating atoms as the impurity.

- `IAO_small/`: Localized orthogonal orbitals based DMET calculations using the 
IAO +PAO partitioning scheme (LO(IAO+PAO)-DMET) with the central lanthanide ion as the impurity.

- `IAO_big/`: LO(IAO+PAO)-DMET calculations with the lanthanide ion and its nearest coordinating atoms as the impurity.

- `IP_small/`: Localized orthogonal orbitals based DMET calculations using the impurity-preserved partition scheme (LO(IP)-DMET) with the central lanthanide ion as the impurity.

- `IP_big/`: LO(IP)-DMET calculations with the lanthanide ion and its nearest
  coordinating atoms as the impurity.

- `Low_small/`: Localized orthogonal orbitals based DMET calculations using the Löwdin partitioning scheme (LO(Löwdin)-DMET) with the central lanthanide ion as the impurity.

- `Low_big/`: LO(Löwdin)-DMET calculations with the lanthanide ion and its nearest coordinating atoms as the impurity.

Individual calculation directories typically contain the Python input script
(`*.py`), PySCF checkpoint files (`*.chk`), job script (`job_submit.sh`) and scheduler output files (`*.out`)
needed to inspect or reproduce the reported results.

## Order of Systems in the Manuscript

The systems appear in the manuscript in the following order:

1. `1Dy`: $[\mathrm{Dy}(\mathrm{Cp}^{\mathrm{ttt}})_2]^+$

2. `2Dy`: $[\mathrm{Dy}\{\mathrm{N}[\mathrm{Si}({}^{i}\mathrm{Pr})_3][\mathrm{Si}({}^{i}\mathrm{Pr})_2\mathrm{C}(\mathrm{CH}_3)=\mathrm{CHCH}_3]\}]^+$

3. `3Ce`: $\mathrm{Ce\text{-}Bp}^{\mathrm{Me}}$

4. `4Ce`: $\mathrm{Ce\text{-}O}_2\mathrm{ip}^{\mathrm{tBu}}$

## Software Notes

The calculation scripts were written for Python workflows based on PySCF and the
in-house `embed_sim` package. To rerun a calculation, install the dependencies
required by the corresponding script, make `embed_sim` importable, For DMET
calculations, copy the corresponding all-electron checkpoint files (`*.chk`) from
the relevant `AE/` directory into the current calculation directory, and then
execute the local `job_submit.sh` or Python input file from that directory.
