---
title: 'Luxum: A parallel C++ library for FDTD solutions to Maxwell equations'
tags:
  - C++
  - parallel
  - finite difference
  - time domain
  - Maxwell
authors:
  - name: Rishabh Sahu
    orcid: 0000-0002-8888-2698
    affiliation: "1"
    email: srishabh@iitk.ac.in
  - name: Mahendra Kumar Verma
    orcid: 0000-0002-3380-4561
    affiliation: "1"
    email: mkv@iitk.ac.in
affiliations:
 - name: Department of Physics, Indian Institute of Technology Kanpur, Kanpur, UP 208016, India
   index: 1
date: 3 August 2018
bibliography: paper.bib
---

# Summary

In our modern world, electromagnetic fields (EM fields) are being employed in many new technologies ranging 
from common devices like cell phones and microwaves to laboratory equipment like lasers. The 
research in further development of these applications and, also, development of new applications 
based on EM fields
---
title: 'Luxum: A parallel C++ library for FDTD solutions to Maxwell equations'
tags:
  - C++
  - parallel
  - finite difference
  - time domain
  - Maxwell
authors:
  - name: Rishabh Sahu
    orcid: 0000-0002-8888-2698
    affiliation: "1"
    email: srishabh@iitk.ac.in
  - name: Mahendra Kumar Verma
    orcid: 0000-0002-3380-4561
    affiliation: "1"
    email: mkv@iitk.ac.in
affiliations:
 - name: Department of Physics, Indian Institute of Technology Kanpur, Kanpur, UP 208016, India
   index: 1
date: 3 August 2018
bibliography: paper.bib
---

# Summary

In our modern world, electromagnetic fields (EM fields) are being employed in many new technologies ranging
from common devices like cell phones and microwaves to laboratory equipment like lasers. The research in further development of these applications and in development of new applications based on EM fields can benefit vastly by studying how light interacts with matter. This allows researchers to pin-point the factors causing errors or inefficiencies. 

One way to study EM waves' interaction with complex-shaped objects is to solve Maxwell's equations numerically. Various attempts have been made to achieve this, solving in both frequency and time domain. After a few decades, these attempts culminated to an algorithm given by Yee in 1966 [@Yee]. The algorithm solves Maxwell's equations in time domain and uses finite difference scheme to evaluate the differentials. 

Since the time Finite Difference Time Domain (FDTD) solutions were introduced, a lot of research has made them complete. FDTD allows a number of features from plane wavefronts sources infinite in one or more direction to materials with complex refractive index. In 1994, Berenger gave an efficient technique to implement effective absorbing boundary condition called Perfectly Matched Layer (PML)[@Berenger]. To top it all, simplicity of FDTD algorithm allows it to be parallelised efficiently. All these factors have made FDTD solutions to Maxwell's equations an excellent choice to study EM field's interaction matter. 
