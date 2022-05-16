function [lamb] = computeElectronWavelength(E0)
%COMPUTEELECTRONWAVELENGTH Compute de Broglie wavelength (in Angstroms)
% of relativistic electron
%   E0 - energy in eV

m = 9.109383*10^-31;
e = 1.602177*10^-19;
c =  299792458;
h = 6.62607*10^-34;
lamb = h/sqrt(2*m*e*E0) ...
    /sqrt(1 + e*E0/2/m/c^2) * 10^10; % wavelength in A
end

