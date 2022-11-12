# polIMParton
Polarized parton distribution functions of the nucleons

polIMParton v1.0

-- pol = polarized

-- IMParton = I'm Parton = Parton from IMP


# Description
This package gives polarized parton distribution functions (pol. PDFs) of the proton
starting from low Q^{2} ~ 0.07 GeV^2, which are based on the analysis to
deep inelastic scattering data applying DGLAP equations with nonlinear
corrections. 

Refs:

"An analysis of polarized parton distribution functions with nonlinear QCD evolution equations", 
C. Han, G. Xie, R. Wang, X. Chen, Nuclear Physics B 985 (2022) 116012;

"Dynamical parton distributions from DGLAP equations with nonlinear corrections", 
R. Wang, X. Chen, Chinese Physics C 41 (2017) 05313.
 


# Usage
The library consists of a C++ class named polIMParton (read ./polIMParton.h
and ./polIMParton.cpp for details). polIMParton has the methods
polIMParton::getPolPDF(Iparton, X, Q2) and polIMParton::getPolPDFError(Iparton, X, Q2),
which are suggested to be called in users' programs. 
Iparton set as -4, -3, -2, -1, 0, 1, 2, 3, 4 corresponds to getting 
cbar, sbar, dbar, ubar, gluon, u, d, s, c
quark/gluon helicity distributions respectively.
The other important method of polIMParton is setDataSet(int).
setDataSet(1) corresponds to use the data set A, 
which is the only data set provided so far.


IMPORTANT!!

./test.cpp gives an example to get polarized up valence quark distribution
in proton and polarized gluon distribution in neutron at high Q^2.

./test.cpp can be modified as users' wants.

**To run the example:**

> cd polIMParton

> cmake .

> make

> ./polPDFtest

