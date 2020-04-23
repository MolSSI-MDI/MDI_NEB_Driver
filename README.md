# MDI Nudged Elastic Band Driver

This repository contains an example driver built using the MolSSI Driver Interface (MDI). It implements a Nudged Elastic Band (NEB) calculation using engines running independent versions of LAMMPS [1]. We modeled our NEB Driver after the LAMMPS neb command, to provide a baseline for testing. This driver uses `n` independent replicas to perform an NEB operation. For each replica, and instance of LAMMPS is run to perform geometry optimizations. Each replica (instance of LAMMPS) requires a LAMMPS input and data file, distributed evenly accross the transition path to correctly operate. Provided with the Driver are a set of LAMMPS files to run the included test script.


## Instructions for Use
This Driver requires LAMMPS, installed through the [MDI fork](https://molssi.github.io/MDI_Library/html/mdi_ecosystem.html#ecosystem_lammps) and the [MDI Driver](https://molssi.github.io/MDI_Library/html/library_page.html). With both codes installed and this Driver downloaded and compiled, there are only a few steps needed to run the tests provided with the Driver.

First, edit the `LAMMPS` and `MDI_NEB_DRIVER` files located within the `NEB_Driver/tests/locations` folder. Edit the file paths to point to the executable of LAMMPS and the compiled driver binary.
Then, simply run the `testing.sh` scrpt in `NEB_Driver/tests/NEB`.

The script has the following command line interface:

`testing.sh num_engines spring_force_constant energy_threshold force_threshold`

  * `num_engines`: the number of engines (replicas) to be used with the driver.
  * `spring_force_constant`: the constant used when determining the spring forces. (see [2] and [3] for a detailed description of the algorithm.)
  * `energy_threshold`: the threshold to check if the energies have converged.
  * `force_threshold`: the threshold to check if the forces have converged.
  










## References
1. S. Plimpton,Fast parallel algorithms for short-range molecular dynam-ics,  Sandia  national  labs.,  albuquerque,  nm  (united  states)  technicalreport, 1993.
2. G. Henkelman, B. P. Uberuaga and H. J ́onsson,The Journal of ChemicalPhysics, 2000,113, 9901–9904.
3. G. Henkelman and H. J ́onsson,The Journal of Chemical Physics, 2000,113, 9978–9985
