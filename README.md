# Multi-Phase-Field-for-FRC-using-Puck-theory

Abaqus implementation of a multi-phase-field fracture model for fiber-reinforced composites using Puck failure theory.

This repository contains the codes and input files used in the paper:

> Pavan Kumar Asur Vijaya Kumar, Rafael Fleischhacker, Aamir Dean, Raimund Rolfes, Heinz E. Pettermann,
> “Revisiting multi-phase field model for FRCs using Puck theory.”
> *Composite Structures*, 372, Article 119549, 2025.
> DOI: https://doi.org/10.1016/j.compstruct.2025.119549

The implementation supports both AT1 and AT2 phase-field fracture models.

---

## Main Features

* Multi-phase-field formulation for fiber-reinforced composites.
* Puck failure theory for distinguishing fiber-dominated and inter-fiber-dominated failure.
* Separate phase-field variables for different failure mechanisms.
* Support for AT1 and AT2 phase-field models.
* Abaqus implementation using user subroutines.
* Example input files corresponding to the numerical examples in the associated publication.
* Material properties and model parameters provided in `pfall.f`.

---

## Software Requirements

The code has been tested with:

* Abaqus 2020
* Abaqus 2022
* Abaqus 2024

A compatible Fortran compiler is required for compiling the Abaqus user subroutines.

---

## Repository Structure

Typical repository contents include:

* `pfall.f`
  Main user subroutine file containing the material properties and model implementation.

* Abaqus input files: `*.inp`
  Input files corresponding to the examples presented in the paper.

* Example folders
  Benchmark cases and numerical examples used to demonstrate the model.

---

## Important Setup Notes

Before running large meshes, check the module `Kvisual`.

The third line contains a fixed-size array:

```fortran
UserVar(70000,16,4)
```

If your mesh contains more than `70000` elements, increase the first dimension to a value larger than the total number of elements in your model.

For example:

```fortran
UserVar(120000,16,4)
```

The selected value should safely exceed the maximum number of elements used in the simulation.

---

## Model Parameters

All required material properties and model parameters are defined in:

```text
pfall.f
```

Before running a simulation, verify:

* elastic properties,
* strength parameters,
* fracture energies,
* phase-field length scales,
* AT1 or AT2 model selection,
* Puck failure parameters,
* mesh size and element definitions,
* state-variable allocation.

It is recommended to start from the provided example input files and modify the parameters gradually.

---

## Running an Example

A typical Abaqus command is:

```bash
abaqus job=<input_file_name> user=pfall.f interactive
```

For example:

```bash
abaqus job=example user=pfall.f interactive
```

Replace `example` with the actual name of the Abaqus input file.

---

## Citation

If you use this code for academic, research, or industrial purposes, please cite the relevant papers.

### Main paper for this repository

```bibtex
@article{AsurVijayaKumar2025RevisitingPuck,
  title   = {Revisiting multi-phase field model for FRCs using Puck theory},
  author  = {Asur Vijaya Kumar, Pavan Kumar and Fleischhacker, Rafael and Dean, Aamir and Rolfes, Raimund and Pettermann, Heinz E.},
  journal = {Composite Structures},
  volume  = {372},
  pages   = {119549},
  year    = {2025},
  doi     = {10.1016/j.compstruct.2025.119549}
}
```

### Related foundational papers

```bibtex
@article{Dean2020PuckPhaseField,
  title   = {A multi phase-field fracture model for long fiber reinforced composites based on the Puck theory of failure},
  author  = {Dean, Aamir and Asur Vijaya Kumar, Pavan Kumar and Reinoso, Jos{\'e} and Gerendt, Carolin and Mahdi, Elsadig and Paggi, Marco and Rolfes, Raimund},
  journal = {Composite Structures},
  volume  = {251},
  pages   = {112446},
  year    = {2020},
  doi     = {10.1016/j.compstruct.2020.112446}
}
```

```bibtex
@article{AsurVijayaKumar2021DelaminationMigration,
  title   = {A multi phase-field-cohesive zone model for laminated composites: Application to delamination migration},
  author  = {Asur Vijaya Kumar, Pavan Kumar and Dean, Aamir and Reinoso, Jos{\'e} and Paggi, Marco},
  journal = {Composite Structures},
  volume  = {276},
  pages   = {114471},
  year    = {2021},
  doi     = {10.1016/j.compstruct.2021.114471}
}
```

---

## Authors and Contact

* Pavan Kumar Asur Vijaya Kumar — GitHub: `Pavan-asur`
  Email: [pavan.kumar@tuwien.ac.at](mailto:pavan.kumar@tuwien.ac.at)
  Email: [asurpavankumar@gmail.com](mailto:asurpavankumar@gmail.com)

* Aamir Dean — GitHub: `aamir-dean`
  Email: [a.dean@sustech.edu](mailto:a.dean@sustech.edu)
  Email: [dr.aamir.dean@gmail.com](mailto:dr.aamir.dean@gmail.com)

* Rafael Fleischhacker

* Raimund Rolfes

* Heinz E. Pettermann

---

## Disclaimer

This repository is research code. The authors provide it to support reproducibility, benchmarking, and further development. Users should carefully verify all model parameters, mesh settings, material definitions, and convergence behavior before using the code for new applications.
