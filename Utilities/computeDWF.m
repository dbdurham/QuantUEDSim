function [B,DWFInt,DWFAmp] = computeDWF(urms,udim,Gmag)
%COMPUTEDWF Compute the isotropic Debye-Waller Factor 
% for diffraction intensities and amplitudes at
% specified reciprocal lattice points
%   urms - RMS atomic displacements (Angstroms)
%   udim - 1, 2, or 3 dimensional atomic displacements
%   Gmag - Reciprocal lattice distances 

B = 8*pi^2*urms^2/udim;
DWFInt = exp(-B*Gmag.^2/4);
DWFAmp = exp(-B*Gmag.^2/2);

end

