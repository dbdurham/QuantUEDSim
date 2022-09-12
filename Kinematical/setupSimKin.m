function sDiff = setupSimKin(options)
%SETUPSIMKIN Set up parameters for kinematical calculations
%   options (optional) - struct containing input variables, listed below

% Default options
E0 = 750e3; % Electron kineatic energy (eV)
uRMS = 0; % 1D rms displacement perpendicular to Bragg planes (Angstroms)
cellMult = 2;
downSampFacCell = 2;
imageSizeCell = 32;

% Overwrite defaults with input options (if provided)
if nargin > 0
    fieldNames = fieldnames(options);
    nFields = length(fieldNames);
    for iField = 1:nFields
        fname = fieldNames{iField};
        [~] = evalc([fname ' = ' num2str(options.(fname))]);
    end
end

% Projected diffraction pattern sampling parameters
downSampFac = cellMult*downSampFacCell;

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
lRange = [-1 1];

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

simType = 'Kinematical';

% Store variables
varsToStore = {'uRMS','E0','lambElec','gamma',...
    'atoms','cellDim','Z','lattice','Vcryst',...
    'hkl','Ghkl','Gmag','dhkl','Gvec',...
    'hRange','kRange','lRange',...
    'Fhkl',...
    'cellMult','imageSizeCell','imageSize',...
    'downSampFac','qxa','qya',...
    'storeSize','qxaStore','qyaStore',...
    'simType'};
nVars = numel(varsToStore);
for iVar = 1:nVars
    [~] = evalc(['sDiff.' varsToStore{iVar} ' = ' varsToStore{iVar}]);
end

end

