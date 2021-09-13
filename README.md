# TransmissionBIE-OpNorms
Plot norms of boundary integral operators associated to Helmholtz transmission problems on circles and spheres

The Matlab function   TransmissionBIE_OpNorms.m   is used to generate the plots in the manuscript
"Spurious quasi-resonances in boundary integral equations for the Helmholtz transmission problem" (2021 - R. Hiptmair, A. Moiola, E.A. Spence).

We solve Helmholtz transmission problems on the unit disc and the unit sphere using Fourier diagonalisation (circular and spherical harmonics).
The eigenvalues of the BIOs are computed in terms of (spherical) Bessel and Hankel functions.
We plot (against the wavenumber k) the norms of the solution operator and the inverses of the 1st/2nd-kind classical/augmented BIOs.

To run the code launch the command 
*  TransmissionBIE_OpNorms

or
*  TransmissionBIE_OpNorms(Dim,ni,no)

where
*  Dim 	is the space dimension, either 2 or 3, 
*  ni	is the PDE coefficient in the bounded domain (1 or 3 in example in the manuscript),  
*  no is the PDE coefficient in the unbounded domain (3 or 1 in example in the manuscript).
  
The figures in the manuscript correspond to 
*  TransmissionBIE_OpNorms(2,3,1)
*  TransmissionBIE_OpNorms(3,3,1)
*  TransmissionBIE_OpNorms(2,1,3)
*  TransmissionBIE_OpNorms(3,1,3)
  
See the comments in the file for details on the formulas used. 
  
Tested on Matlab R2019b.

