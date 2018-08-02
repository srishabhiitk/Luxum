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
   * Point Sources
   * Plane Wavefront Sources
 * Boundary Conditions:
   * Perfect Electric Conductor (PEC) boundary which relfects any scattered back into the domain
   * Absorbing Boundary Condition using Perfectly Matched Layer

 Further work on development of this library is going on. We welcome collaborations in developing the code further. Interested people can directly write to [Prof. Mahendra Verma](mailto:mkv@iitk.ac.in).

 Code is being developed by [Mr. Rishabh Sahu](https://scholar.google.co.in/citations?user=iEj0p54AAAAJ&hl=en).


# Example

Below are examples of different features mentioned above and some other applications. These examples show the usefulness of the developed library in certain applications like simulation of photonic crystals. As the code is parallel, it can be easily used to simulate huge integrated optics chips.

---

## Two dimensional simulation

FDTD simulation in 2D where a pulse with plane wavefront hits a square crystal of refractive index 3.3.

 [![Alt text for your video](https://img.youtube.com/vi/eJ4aiyMKMZY/0.jpg)](http://www.youtube.com/watch?v=eJ4aiyMKMZY)

 ---

## Three dimensional simulation

FDTD simulation in 3D where a pulse with plane wavefront hits a cubic crystal of refractive index 3.3.

[![Alt text for your video](https://img.youtube.com/vi/tloQzR-SyAk/0.jpg)](http://www.youtube.com/watch?v=tloQzR-SyAk)

---

## Perfectly Matched Layer

Below simulation shows a demostration of working of PML. It is a FDTD simulation of a radiating source. On the left side, there is no PML implementation. Hence, the light reflects back and superposes with the existing light. On the right, there is a PML on LHS and, therefore, it does not reflect back.

[![Alt text for your video](https://img.youtube.com/vi/-eT4g0tGn50/0.jpg)](http://www.youtube.com/watch?v=-eT4g0tGn50)

The effect of PML  can also be shown in three dimensions. Below is the simulation which shows effect of PML when a plane wavefront interacts with a cube of relative electric permeability of 6.

[![Alt text for your video](https://img.youtube.com/vi/hs9v3FfOKYk/0.jpg)](http://www.youtube.com/watch?v=hs9v3FfOKYk)


---

## Lens

FDTD simulation in 2D where a pulse with plane wavefront hits a lens of refractive index 1.5.

[![Alt text for your video](https://img.youtube.com/vi/1KqlazgbGBs/0.jpg)](http://www.youtube.com/watch?v=1KqlazgbGBs)

Light can also be made incident with a oblique angle of incidence.

[![Alt text for your video](https://img.youtube.com/vi/1rm_AgInxyE/0.jpg)](http://www.youtube.com/watch?v=1rm_AgInxyE)


---

## Distorted Lens

FDTD simulation in 2D where a pulse with plane wavefront hits a distorted lens of refractive index 1.5.

[![Alt text for your video](https://img.youtube.com/vi/SIFs8SABcWg/0.jpg)](http://www.youtube.com/watch?v=SIFs8SABcWg)


---

## Photonic Crystal Waveguide

FDTD simulation of a photonic crystal waveguide. It is made using Silicon square blocks arranged in a periodic lattice. The relative electric permittivity of silicon is taken to be 12. The gap is filled with air or vacuum.
The video shows confinement of light in a waveguide.

[![Alt text for your video](https://img.youtube.com/vi/2S5wH669rdw/0.jpg)](http://www.youtube.com/watch?v=2S5wH669rdw)

Coupling of waveguides can also be observed if another waveguide very close to existing waveguide is put.

[![Alt text for your video](https://img.youtube.com/vi/ahHDnKWI5-k/0.jpg)](http://www.youtube.com/watch?v=ahHDnKWI5-k)


The coupling changes if the distance between two waveguides is changed.

[![Alt text for your video](https://img.youtube.com/vi/HS1VBHk3l1o/0.jpg)](http://www.youtube.com/watch?v=HS1VBHk3l1o)


---

## Complex 3D Object Simulation

FDTD simulation of a plane wavefront hitting a metallic airplane. The boundary is absorbing using Perfectly Matched Layer.

[![Alt text for your video](https://img.youtube.com/vi/bildTtxlKsY/0.jpg)](http://www.youtube.com/watch?v=bildTtxlKsY)
