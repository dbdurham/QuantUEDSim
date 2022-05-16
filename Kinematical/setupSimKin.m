function sDiff = setupSimKin()
%SETUPSIMKIN Set up parameter struct for kinematical diffraction
%simulations
%   

% Input variables
E0 = 750e3; %eV
uRMS = 0.0894; % 1D rms displacement perpendicular to Bragg planes (Angstroms)

lambElec = computeElectronWavelength(E0); % Angstroms
gamma = computeLorentzFactor(E0); % 

%% Generate crystal structure and reciprocal lattice

% Generate crystal structure: unit cell and atomic coordinates
options.cellMult = 1;
[atoms,cellDim,lattice,uvwInit] = wyckoffGold(options);
Z = atoms(:,4); % Atomic number
nAtoms = size(lattice,1); % Number of atoms in unit cell
Vcryst = prod(cellDim); % Unit cell volume (Angstroms^3)

% Generate mesh of hkl
hRange = [-8 8];
kRange = [-8 8];
lRange = [-1 1];

[hkl,Ghkl,Gmag,dhkl] = generateReciprocalLattice(...
    uvwInit,hRange,kRange,lRange);

%% Calculate structure factors
% NOTE: scattering factors assumed independent of beam energy... an
% approximation often made for relativistic electrons

Fhkl = computeStructureFactors(lattice,Z,hkl,Gmag,E0,...
    'Moliere',false);

%% Store variables
% Store variables
varsToStore = {'uRMS','E0','lambElec','gamma',...
    'atoms','cellDim','Z','lattice','Vcryst',...
    'hkl','Ghkl','Gmag','dhkl',...
    'hRange','kRange','lRange',...
    'Fhkl'};
nVars = numel(varsToStore);
for iVar = 1:nVars
    [~] = evalc(['sDiff.' varsToStore{iVar} ' = ' varsToStore{iVar}]);
end

end

