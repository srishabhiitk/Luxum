**Luxum** is a parallel, C++ based Finite Difference Time Domain solver for
Maxwell's equation.


# Introduction
The interest in time-stepping solutions of electromagnetic fields using Maxwell’s equations was spurred by development of microwave radar technology and it remained the prime motivation until 1990. During this time, development of stealth aircrafts and even better radar technology was of paramount importance. However, as time passed many other applications came to light. They ranged from designing photonic crystals and low power optical switches to early stage detection of breast cancer and observing cell phone radiation interaction with human head.

FDTD solutions have been around for several decades and various softwares from free open-source licenses to expensive ones exist now. These softwares cover a variety of range of applications, however, an FDTD library that is designed around parallelizability to tackle extremely large simulations using supercomputers is still missing.

In this project, we develop a C++ library that feature near ideal parallelizability along with other general FDTD features such as support for plane-wave sources and PML absorbing boundary condition. We have also implemented support for parallel file reading and writing for efficient IO using HDF5 library. Further, many companions codes are developed that work along with library for plotting, generating animations and reading images to input material shapes.

For further details on code developement of this project, please refer to [FDTD Solver library](http://http://turbulencehub.org/index.php/codes/fdtd-solver/).


# Code
This project is dedicated to develop an object-oriented parallel C++ code that can simulate electromagnetic waves using Finite Difference Time Domain scheme. At the current stage, the developed library supports a number of features,
 * Two and three dimensional simulations
   (with spatial material characterization)
 * Sources:
 * * Point Sources
 * * Plane Wavefront Sources
 * Boundary Conditions:
 * * Perfect Electric Conductor (PEC) boundary which relfects any scattered back into the domain
 * * Absorbing Boundary Condition using Perfectly Matched Layer

 Further work on development of this library is going on. We welcome collaborations in developing the code further. Interested people can directly write to [Prof. Mahendra Verma](mailto:mkv@iitk.ac.in).

 Code is being developed by [Mr. Rishabh Sahu](https://scholar.google.co.in/citations?user=iEj0p54AAAAJ&hl=en).
