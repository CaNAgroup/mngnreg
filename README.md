# mngnreg
minimal-norm solution of a nonlinear least-squares problem.

Package attached to the paper:

_Federica Pes and Giuseppe Rodriguez_,
**The MNGNREG toolbox for the regularized solution of nonlinear least-squares problems**,
Submitted, 2025

---

Authors
-------

- Federica Pes, Giuseppe Rodriguez,
       University of Cagliari, Italy

> Email: federica.pes@unica.it, rodriguez@unica.it

---

Description
-----------

MNGNREG: minimal-norm solution of a nonlinear least-squares problem.

Release: 2.0, April 2025.

Programming language: Matlab 24.1 (R2024a).

License: open source.

This Matlab package implements the algorithms described in the papers:

F. Pes and G. Rodriguez.
The minimal-norm Gauss-Newton method and some of its regularized variants.
Electron. Trans. Numer. Anal., 53:459-480, 2020.

F. Pes and G. Rodriguez.
A doubly relaxed minimal-norm Gauss-Newton method for underdetermined nonlinear least-squares problems.
Appl. Numer. Math., 171:233-248, 2022.
  
F. Pes.
Truncated Minimal-Norm Gauss-Newton Method applied to the inversion of FDEM data.
Computational Science and Its Applications - ICCSA 2023 Workshops,
Lecture Notes in Computer Science, vol 14108, 641-658, 2023.

---

Installation
------------

After downloading the package, to use it, the user can add the routines directly into the working directory or can add mngnreg to the Matlab search path by the "addpath" command.

---

PACKAGE USE
-----------

Some simple drivers are provided.

---

Package structure
-----------------

The package contains the following files:
* README.md   : this file
* tmngn       : minimal-norm Gauss-Newton regularized by TSVD
* tmlngn      : minimal-norm Gauss-Newton regularized by TGSVD
* tikgn       : minimal-norm Gauss-Newton regularized by Tikhonov in standard form
* tiklgn      : minimal-norm Gauss-Newton regularized by Tikhonov in general form
* driverUnder : test program for tmngn.m with beta=1
* driverOver  : test program for tmlngn.m with beta=1
* driverMNGN2 : test program for tmngn.m
* funUnder    : generates a simple underdetermined nonlinear problem
* funOver     : generates a simple overdetermined nonlinear problem
* nonlinfun1  : generates an underdetermined nonlinear problem
* get_l       : compute discrete derivative operators (from Regutools)
* jack        : approximate Jacobian matrix by finite differences (from FDEMtools)

This package uses the routines:

get_l.m from the Regularization Tools toolbox by P.C. Hansen 
https://www.imm.dtu.dk/~pcha/Regutools/

jack.m from the FDEMtools by G.P. Deidda, P. Díaz De Alba, C. Fenu, G. Lovicu, and G. Rodriguez
https://bugs.unica.it/cana/software/#FDEMtools3

---

Package updates
---------------

The authors may indicate a web page where package updates will be available
[link](https://github.com/CaNAgroup/mngnreg)
