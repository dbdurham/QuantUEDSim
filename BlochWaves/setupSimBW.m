function sDiff = setupSimBW(options)
%SETUPSIMBW Set up the Bloch Wave calculation
%   options (optional) - struct containing input variables, listed below

% Default options
E0 = 750e3; % eV beam kinetic energy    
uRMS = 0.0894; % 1D rms displacement across Bragg planes (Angstroms)
GxyThresh = 4.5; % in-plane reciprocal space threshold (inv Angstroms)
sThresh = 0.1; % Excitation error threshold (inv Angstroms)
cellMult = 2;
downSampFacCell = 2;
imageSizeCell = 32;
UgThresh = 1e-5;

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
sigma = computeInteractionParameter(E0); % rad/(V*Angstroms)


%% 1. Generate crystal structure and reciprocal lattice

% Generate crystal structure: unit cell and atomic coordinates
options.cellMult = 1;
[atoms,cellDim,lattice,uvwInit] = wyckoffGold(options);
Z = atoms(:,4); % Atomic number
Vcryst = prod(cellDim); % unit cell volume (Angstroms^3)

% Generate mesh of hkl
hRange = [-50 50];
kRange = [-50 50];
lRange = [-30 30]; 
[hkl,Ghkl,Gmag,dhkl,Gvec] = generateReciprocalLattice(...
    uvwInit,hRange,kRange,lRange);

hLen = diff(hRange)+1;
kLen = diff(kRange)+1;
lLen = diff(lRange)+1;

% Diffraction pattern sampling
imageSize = imageSizeCell.*ones(1,2); 
storeSize = imageSize./downSampFac;
pixelSize = cellMult*cellDim(1:2)./imageSize;
[qxa,qya] = makeFourierCoords(imageSize,pixelSize);
[qxaStore,qyaStore,storeMask] = downsampleFourierCoords(...
    qxa,qya,downSampFac);

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

simType = 'Bloch Waves';

% Store variables
varsToStore = {'uRMS','E0','lambElec','sigma',...
    'atoms','cellDim','Z','lattice','Vcryst',...
    'hkl','Ghkl','Gmag','Gvec','dhkl',...
    'hRange','kRange','lRange','hLen','kLen','lLen',...
    'U_G_0K','U_G','GxyThresh','sThresh','UgThresh',...
    'cellMult','imageSizeCell','imageSize',...
    'downSampFac','qxa','qya',...
    'storeSize','qxaStore','qyaStore',...
    'simType'};
nVars = numel(varsToStore);
for iVar = 1:nVars
    [~] = evalc(['sDiff.' varsToStore{iVar} ' = ' varsToStore{iVar}]);
end

end

