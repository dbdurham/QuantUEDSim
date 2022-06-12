function sDiff = setupSimKin()
%SETUPSIMKIN Set up parameter struct for kinematical diffraction
%simulations
%   

% Input variables
E0 = 750e3; %eV
uRMS = 0.0894; % 1D rms displacement perpendicular to Bragg planes (Angstroms)

% Projected diffraction pattern sampling parameters
cellMult = 2;
imageSizeCell = 64;
downSampFac = cellMult;

% Beam physical constants
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
lRange = [-5 5];

[hkl,Ghkl,Gmag,dhkl,Gvec] = generateReciprocalLattice(...
    uvwInit,hRange,kRange,lRange);

% Diffraction pattern sampling
imageSize = imageSizeCell.*ones(1,2); 
storeSize = imageSize./downSampFac;
pixelSize = cellMult*cellDim(1:2)./imageSize;
[qxa,qya] = makeFourierCoords(imageSize,pixelSize);
% Downsampling
[qxaStore,qyaStore,storeMask] = downsampleFourierCoords(...
    qxa,qya,downSampFac);

%% Calculate structure factors
Fhkl = computeStructureFactors(lattice,Z,hkl,Gmag,E0,...
    'Moliere',false);

%% Store variables
% Store variables
varsToStore = {'uRMS','E0','lambElec','gamma',...
    'atoms','cellDim','Z','lattice','Vcryst',...
    'hkl','Ghkl','Gmag','dhkl','Gvec',...
    'hRange','kRange','lRange',...
    'Fhkl',...
    'cellMult','imageSizeCell','imageSize',...
    'downSampFac','qxa','qya',...
    'storeSize','qxaStore','qyaStore'};
nVars = numel(varsToStore);
for iVar = 1:nVars
    [~] = evalc(['sDiff.' varsToStore{iVar} ' = ' varsToStore{iVar}]);
end

end

