function [sigma] = computeInteractionParameter(E0)
%COMPUTEINTERACTIONPARAMETER Computes the electron interaction parameter in
%rad/(V*Angstroms) for electrons of given kinetic energy, including
%relativistic correction
%   E0 = kinetic energy (eV)

lambElec = computeElectronWavelength(E0);
m = 9.11e-31; % Electron mass (kg)
c = 3e8; % Speed of light (m/s)
e = 1.602e-19; % Charge of electron (C)

sigma = (2*pi/lambElec/E0) ...
    *(m*c^2+e*E0)/(2*m*c^2+e*E0); % Interaction parameter (rad / (V*A))

end

