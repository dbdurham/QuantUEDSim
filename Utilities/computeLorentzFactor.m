function [gamma] = computeLorentzFactor(E0)
%COMPUTELORENTZFACTOR Compute electron Lorentz factor for given kinetic
%energy
%   E0 - energy in eV

m = 9.109383*10^-31;
e = 1.602177*10^-19;
c =  299792458;
gamma = 1+(e*E0/m/c^2);

end

