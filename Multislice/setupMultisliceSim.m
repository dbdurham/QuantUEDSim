%% Setup dynamical diffraction simulation for 
function [sDiff] = setupMultisliceSim(options,wyckOptions)
    
showKin = true; % calculate and show "kinematical pattern" from one UC

% Wyckoff options
sDiff.cellMult = 1;

% Incident beam angles
sDiff.theta_x = 0;
sDiff.theta_y = 0;

% pixel size options: either specify desired size, or the image size
% directly
% sDiff.pixelSize = 0.02;           % approx realspace pixel size (Angstroms).
% f = 4;
% sDiff.imageSize = round(sDiff.cellDim(1:2)/sDiff.pixelSize./f).*f;
sDiff.imageSizeCell = 256;

sDiff.potBound = 3;       % Radial distance to integrate atomic potentials. (real space, Angstroms)
% Note that a higher potBound is needed for DW damping which broadens the
% pots! Check for convergence at each temperature!
% sDiff.fullImagePot = true;
sDiff.E0 = 750e3;  % Microscope accel voltage in volts

% Thermal motion
sDiff.uRMS = 0.0894; % (1D) rms displacement perpendicular to diff planes (Angstroms), 0.165 at 300 K 
sDiff.numFP = 1; %(for 0 K sims or DW damping, FP = 1); 
sDiff.DWDamping = true; % apply DW damping to the crystal potentials... do not use with FPs!
sDiff.computePot = false; % Compute a new potential for each UC
% (set to false to use frozen phonons instead)
sDiff.anisoCryst = 0;

% Manual definition of slice thicknesses
% sRefine.numSlices = 2;  % number of slices along beam direction
% sRefine.sliceThickness = [0.5 0.5]; % relative size of slices (optional
% for slices of equal thickness)
useAutoSlice = true;
minThickness = 0; 

sDiff.flagUseAntiAliasing = true;
sDiff.AntiAliasingMode = '1/2'; % '1/2' band limits propagator to kmax/2, 
% '2/3' band limits prop to (2/3)*kmax and expPot to (2/3)*kmax

% Output size
sDiff.qRange = 4; % q range to store in *output* (A^-1)

% Overwrite defaults with input options (if provided)
if nargin > 0
    fieldNames = fieldnames(options);
    nFields = length(fieldNames);
    for iField = 1:nFields
        fname = fieldNames{iField};
        sDiff.(fname) = options.(fname);
    end
end

if sDiff.DWDamping
    sDiff.M = 2*pi^2*sDiff.uRMS^2; % DW damping factor on the potentials -- in future, could be atom-dependent
    sDiff.uFP = 0; % rms displacement of atoms in frozen phonon simulations
else
    sDiff.M = 0;
    sDiff.uFP = sDiff.uRMS; 
end
 
sDiff.downSampFac = sDiff.cellMult; % only keep every downSampFac pixels... Is a downSampFac^2 reduction in memory
% NOTE: this is not equivalent to the q range of the simulation during the
% computation, only what is stored after completion!

wyckOptions.cellMult = sDiff.cellMult;

% Generate simulation cell: parameters and atomic coordinates
[sDiff.atoms,sDiff.cellDim] = wyckoffGold(wyckOptions);

% Compute final pixel size
sDiff.imageSize = sDiff.imageSizeCell*sDiff.cellMult;
sDiff.imageSize = sDiff.imageSize.*ones(1,2); % make 2D if scalar
sDiff.pixelSize = sDiff.cellDim(1:2) ./ sDiff.imageSize;

% slices
if useAutoSlice % auto-slice based on atomic positions
    if ~isfield(sDiff,'numSlices')
        sDiff.numSlices = Inf;
    end
    sDiff.sliceThickness = computeSlices(sDiff.atoms(:,3),...
        sDiff.cellDim(3),minThickness,sDiff.numSlices);
    % update numSlices
    sDiff.numSlices = numel(sDiff.sliceThickness);
    % assign atoms to slices
    tAccum = cumsum(sDiff.sliceThickness);
    sDiff.sliceIndex = zeros(size(sDiff.atoms(:,3)));
    numAtoms = size(sDiff.atoms,1);
    for iAtom = 1:numAtoms
        sDiff.sliceIndex(iAtom) = ...
            find(sDiff.atoms(iAtom,3) < tAccum,1,'first');
    end
elseif isfield(sDiff,'sliceThickness') % pre-defined slices
    sDiff.sliceThickness = sDiff.cellDim(3)...
        .*sDiff.sliceThickness...
        ./sum(sDiff.sliceThickness);
    % assign atoms to slices
    tAccum = cumsum(sDiff.sliceThickness);
    sDiff.sliceIndex = zeros(size(sDiff.atoms(:,3)));
    numAtoms = size(sDiff.atoms,1);
    for iAtom = 1:numAtoms
        sDiff.sliceIndex(iAtom) = ...
            find(sDiff.atoms(iAtom,3) < tAccum,1,'first');
    end
else % equal slices
    sDiff.sliceThickness = sDiff.cellDim(3) / sDiff.numSlices;
    % assign atoms to slices
    sDiff.sliceIndex = min(max(round( ...
        sDiff.atoms(:,3) / sDiff.cellDim(3) * sDiff.numSlices + 0.5),...
        1),sDiff.numSlices);
end

% Calculate wavelength and electron interaction parameter
sDiff.lambda = computeElectronWavelength(sDiff.E0);
sDiff.sigma = computeInteractionParameter(sDiff.E0); % V^-1 Ang^-1

% Coords
qx = makeFourierCoords(sDiff.imageSize(1),sDiff.pixelSize(1));
qy = makeFourierCoords(sDiff.imageSize(2),sDiff.pixelSize(2));
[sDiff.qya,sDiff.qxa] = meshgrid(qy,qx);
q2 = sDiff.qxa.^2 + sDiff.qya.^2;

% Reduced coords for output storage (to save storage space)
if sDiff.qRange > max(sDiff.qxa(:)) || sDiff.qRange > max(sDiff.qya(:))
    sDiff.storeMask = true(size(sDiff.qxa));
else
    sDiff.storeMask = (abs(sDiff.qxa) < sDiff.qRange) & ...
        (abs(sDiff.qya) < sDiff.qRange);
%     xLim = find(fftshift() > sDiff.qRange,1,'first');
%     yLim = find(fftshift(sDiff.qya(1,:)) > sDiff.qRange,1,'first');
%     sDiff.yCrop = (sDiff.imageSize(2)-yLim+1):yLim;
%     sDiff.xCrop = (sDiff.imageSize(1)-xLim+1):xLim;
end
% Downsampling
dqx = sDiff.qxa(2,1)-sDiff.qxa(1,1);
dqy = sDiff.qya(1,2)-sDiff.qya(1,1);
pxa = round(sDiff.qxa/dqx); % pixels away from center
pya = round(sDiff.qya/dqy);
if sDiff.downSampFac > 1
    sDiff.storeMask = sDiff.storeMask & ...
         mod(pxa,sDiff.downSampFac) == 0 & ...
         mod(pya,sDiff.downSampFac) == 0;
end
sDiff.storeLenX = sum(sDiff.storeMask(:,1));
sDiff.storeLenY = sum(sDiff.storeMask(1,:));
sDiff.qxaStore = reshape(sDiff.qxa(sDiff.storeMask),...
    [sDiff.storeLenY,sDiff.storeLenX]);
sDiff.qyaStore = reshape(sDiff.qya(sDiff.storeMask),...
    [sDiff.storeLenY,sDiff.storeLenX]);
% sDiff.qxaStore = fftshift(sDiff.qxa);
% sDiff.qxaStore = fftshift(sDiff.qxaStore(sDiff.xCrop,sDiff.yCrop));
% sDiff.qyaStore = fftshift(sDiff.qya);
% sDiff.qyaStore = fftshift(sDiff.qyaStore(sDiff.xCrop,sDiff.yCrop));

% Propagator
switch sDiff.AntiAliasingMode
    case '1/2'
        sDiff.qMask = double(q2 < (0.5*min(max(abs(qx)),max(abs(qy))))^2);
    case '2/3'
        sDiff.qMask = double(q2 < ((2/3)*min(max(abs(qx)),max(abs(qy))))^2);
end
t = reshape(sDiff.sliceThickness,...
    [1 1 length(sDiff.sliceThickness)]);
if sDiff.flagUseAntiAliasing == false
    sDiff.prop =  exp( ...
        (-1i*pi*sDiff.lambda*t).*q2);
else
    sDiff.prop = sDiff.qMask .* exp( ...
        (-1i*pi*sDiff.lambda*t).*q2);
end
sDiff.propConst = (-1i*pi*sDiff.lambda*t).*q2;

% Input beam
sDiff.EWInit = zeros(sDiff.imageSize);
if isfield(sDiff,'qSigma') & sDiff.qSigma > 0
    sDiff.EWInit(:) = exp(-0.5*(sDiff.qxa.^2 + sDiff.qya.^2)/sDiff.qSigma^2)...
        ./ (sDiff.qSigma*sqrt(2*pi));
else
    sDiff.EWInit(1,1) = 1;
end

% Construct initial potential
tic
sDiff = computePotential(sDiff);
disp('Computed potentials: ')
toc


if showKin

    % Calculate the diffraction pattern from 1 UC (approximately kinematical)  
    sDiff.EWKin = ones(sDiff.imageSize);
    for a0 = 1:sDiff.numSlices
        % Account for uneven slice propagator if needed
        if size(sDiff.prop,3) > 1
            prop = sDiff.prop(:,:,a0);
        else
            prop = sDiff.prop;
        end
        sDiff.EWKin(:) = ifft2(fft2(sDiff.EWKin .* ...
            exp(1i*sDiff.sigma*sDiff.pot(:,:,a0))) .* ...
            prop);
    end
    % diffraction pattern
    sDiff.imageDiffKin = abs(fft2(sDiff.EWKin)).^2;


    % plotting
    Ip = fftshift(sDiff.imageDiffKin).^0.5;
    Isort = sort(Ip(:));
    figure(56)
    clf
    imagesc(Ip)
    axis off
    
  
    colormap(violetFire(256))
%     colormap(simColourMap(2));
    colorbar
    caxis([0 Isort(end-1)])
    
end

end