function sDiff = setupSimBW()
%SETUPSIMBW Set up the Bloch Wave calculation
%   Detailed explanation goes here

% Inputs
E0 = 750e3; %eV
uRMS = 0.0894; % 1D rms displacement across Bragg planes (Angstroms)

lambElec = computeElectronWavelength(E0); % Angstroms
sigma = computeInteractionParameter(E0); % rad/(V*Angstroms)


%% 1. Generate crystal structure and reciprocal lattice

% Generate crystal structure: unit cell and atomic coordinates
options.cellMult = 1;
[atoms,cellDim,lattice,uvwInit] = wyckoffGold(options);
Z = atoms(:,4); % Atomic number
Vcryst = prod(cellDim); % unit cell volume (Angstroms^3)

% Generate mesh of hkl
hRange = [-26 26];
kRange = [-26 26];
lRange = [-2 2]; 
[hkl,Ghkl,Gmag,dhkl] = generateReciprocalLattice(...
    uvwInit,hRange,kRange,lRange);

hLen = diff(hRange)+1;
kLen = diff(kRange)+1;
lLen = diff(lRange)+1;

%% 2. Compute library of 3D specimen potential Fourier components U_G

% Compute structure factors
Fhkl = computeStructureFactors(lattice,Z,hkl,Gmag,E0,...
    'Born',false);

% Scattering potential Fourier components
U_G_0K = (sigma/(pi*lambElec))*(47.86/Vcryst).*Fhkl;
% Debye-Waller factor
[~,~,DWFAmp] = computeDWF(uRMS,1,Gmag);
U_G = U_G_0K.*DWFAmp;
% Force U_G to be real (correct for centrosymmetric crystals)
U_G = real(U_G); 

% Store variables
varsToStore = {'uRMS','E0','lambElec','sigma',...
    'atoms','cellDim','Z','lattice','Vcryst',...
    'hkl','Ghkl','Gmag','dhkl',...
    'hRange','kRange','lRange','hLen','kLen','lLen',...
    'U_G_0K','U_G'};
nVars = numel(varsToStore);
for iVar = 1:nVars
    [~] = evalc(['sDiff.' varsToStore{iVar} ' = ' varsToStore{iVar}]);
end

end

