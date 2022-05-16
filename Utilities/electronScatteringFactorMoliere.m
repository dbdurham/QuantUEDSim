function f = electronScatteringFactorMoliere(Z,q,E0,fparams)
%ELECTRONSCATTERINGFACTORMOLIERE Compute the electron scattering factor
%using the "Moliere approximation," which is the Fourier transform of the
%weak phase object approximation to the real-space scattering process
%(Kirkland Advanced Computing in EM 3rd edition, eq. 5.18.)
%   Z - atomic number
%   q - scattering magnitude
%   E0 - kinetic energy (eV)
%   fparams - scattering factor parameter table from Kirkland

% Relativistic electron parameters
lambElec = computeElectronWavelength(E0); % Angstroms
sigma = computeInteractionParameter(E0); % (rad / (V*A))

% Projected potential
a0 = 0.5292; % Bohr radius in Angstroms
eVA = 14.4; % Volt-Angstroms
% term1 = 2*pi^2*a0*eVA;
% term2 = 2*pi^(5/2)*a0*eVA;
term1 = 4*pi^2*a0*eVA;
term2 = 2*pi^2*a0*eVA;

ap = fparams(Z,:);
potProj = @(r) term1*( ...
    ap(1)*besselk(0,2*pi*sqrt(ap(2))*r) ...
    + ap(3)*besselk(0,2*pi*sqrt(ap(4))*r) ...
    + ap(5)*besselk(0,2*pi*sqrt(ap(6))*r)) ...
    + term2*( ...
    ap(7)/ap(8)*exp(-pi^2/ap(8)*r.^2) ...
    + ap(9)/ap(10)*exp(-pi^2/ap(10)*r.^2) ...
    + ap(11)/ap(12)*exp(-pi^2/ap(12)*r.^2));

% Integrand
intFunc = @(r) besselj(0,2*pi*q.*r).*(1-exp(1i*sigma*potProj(r))).*r;

% Radial integration
f = (2*pi*1i/lambElec) * integral(intFunc,0,Inf);


end

