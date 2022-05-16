function [EW,EWstore,sDiff,IntStore] = calcDiff(sDiff,coefs,expPot)

% Colin Ophus - 2020 April
% Modified by Dan Durham and Khalid Siddiqui
% Calculate multislice diffraction pattern

% sDiff is a struct containing the parameters needed to perform the
% diffraction calculation

% coefficients
numUCs = coefs(1); % can be fractional
% B = max(coefs(2),0); % DW Factor -- NOT USED
theta_x = coefs(3); %  hor tilt
theta_y = coefs(4); % vert tilt

sDiff.theta_x = theta_x;
sDiff.theta_y = theta_y;

% Compute potential at prescribed sample tilt
sDiff = computePotential(sDiff);

% Option to recompute pots for each unit cell
computePot = sDiff.computePot;

% store the exit waves from every unit cell (UC)
flagStoreEWs = true;
% option to store single cell real-space wave functions
if ~isfield(sDiff,'storeRealSpaceCellEWs')
    storeRealSpaceCellEWs = false;
else
    storeRealSpaceCellEWs = sDiff.storeRealSpaceCellEWs;
end
% sample the exit waves at sub-UC scale
if ~isfield(sDiff,'numSubUCs')
    numSubUCs = 1; % number of sub-UCs within each UC
else
    numSubUCs = sDiff.numSubUCs;
end
% set to 1 to only sample the UCs
numSlicesSubUC = ceil(sDiff.numSlices/numSubUCs);
numSlices = sDiff.numSlices;

% Recompute propagator with tilts
% For uneven slices, sliceThickness is a vector with numSlices
% elements containing the thicknesses of each slice
t = reshape(sDiff.sliceThickness,...
        [1 1 length(sDiff.sliceThickness)]);
propType = 'Parabolic';
switch propType
    case 'Parabolic'
        sDiff.prop = exp( ...
            sDiff.propConst ...
            + (2i*pi*t.* ...
            (tan(theta_x)*sDiff.qxa + tan(theta_y)*sDiff.qya)) ...
            );
    case 'Spherical' 
        % Spherical wave propagator
        k = 1/sDiff.lambda;
        phi = atan2(theta_y,theta_x)+pi;
        theta_tot = atan(sqrt(tan(theta_x)^2 + tan(theta_y)^2)); 
        kx = k*cos(phi)*sin(theta_tot);
        ky = k*sin(phi)*sin(theta_tot);
        kz = k*cos(theta_tot);
        % Direct Fourier space definition 
        kMinkPlanar = sqrt((k^2 - ((kx + sDiff.qxa).^2 + (ky+sDiff.qya).^2)));
        Sq = (kz -  kMinkPlanar); % Assume orthogonal crystal system
%         sDiff.prop = (kz ./ kMinkPlanar) ...
%             .* exp(-2i*pi*t.*Sq);
        sDiff.prop  = exp(-2i*pi*t.*Sq);
end
if sDiff.flagUseAntiAliasing
    sDiff.prop = sDiff.qMask .* sDiff.prop;
end

% Determine if slices are empty or not (could do in setup phase)
isEmptySlice = false(sDiff.numSlices,1);
for iSlice = 1:numSlices
    isEmptySlice(iSlice) = ~any(reshape(sDiff.pot(:,:,iSlice),1,[]));
end

if nargin < 3
    expPot = single(sDiff.expPot);
end
prop = single(sDiff.prop);
EW = single(sDiff.EWInit);

% Convert quantities to singles to reduce memory and performance
% requirements
% useSingle = true;
% fieldList = {'prop','pot','qxa','qya','sigma'};
% if useSingle
%     for iField = 1:length(fieldList)
%         sDiff.(fieldList{iField}) = single(sDiff.(fieldList{iField}));
%     end
%     EW = single(EW);
% end

if flagStoreEWs
    if storeRealSpaceCellEWs
        EWstore = zeros([size(EW)./sDiff.cellMult ...
            round(numUCs*numSubUCs)],'single');
    else
        EWstore = zeros(...
            [sDiff.storeLenX,sDiff.storeLenY,round(numUCs*numSubUCs)],'single');  
    end
    IntStore = zeros(1,round(numUCs*numSubUCs),'single');
else
    EWstore = [];
end

for a0 = 1:ceil(numUCs)
    if a0 == ceil(numUCs) % last 
        numSlices = round((numUCs-(a0-1))*sDiff.numSlices);
    end
    if computePot
        sDiff = computePotential(sDiff);
        expPot = single(exp(1i*sigma*sDiff.pot));
    end
    
    for a1 = 1:numSlices
        
        % Interact and propagate
        if ~isEmptySlice(iSlice) % if potential slice is nonzero, full interaction
            EW(:) = prop(:,:,a1) .* fft2(ifft2(EW) .* ...
                    expPot(:,:,a1,randi(sDiff.numFP)));
        else % if potential slice is empty, just propagate
            EW(:) = prop(:,:,a1) .* EW;
        end
        
        % Store output
        if flagStoreEWs && mod(a1,numSlicesSubUC) == 0
            indStore = (a0-1)*numSubUCs+(a1/numSlicesSubUC);
            if storeRealSpaceCellEWs
                EWr = ifft2(EW.*numel(EW)); % Real space EW normalized
                % to have an average value of 1
                EWstore(:,:,indStore) = ...
                    EWr(1:size(EW,1)/sDiff.cellMult,...
                    1:size(EW,2)/sDiff.cellMult);
            else
                EWstore(:,:,indStore) = reshape(EW(sDiff.storeMask),...
                    [sDiff.storeLenX,sDiff.storeLenY]);
            end
        end

    end
     
end


end