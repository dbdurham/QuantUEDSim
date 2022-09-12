function [EW,EWstore,sDiff] = calcDiffMSGPU(sDiff,coefs,expPot)

% Colin Ophus - 2020 April
% Modified by Dan Durham and Khalid Siddiqui
% Calculate multislice diffraction pattern using GPU computing

% sDiff is a struct containing the parameters needed to perform the
% diffraction calculation

computePot = sDiff.computePot;

% store the exit waves from every unit cell (UC)
flagStoreEWs = true;
% sample the exit waves at sub-UC scale
if ~isfield(sDiff,'numSubUCs')
    numSubUCs = 1; % number of sub-UCs within each UC
else
    numSubUCs = sDiff.numSubUCs;
end
% set to 1 to only sample the UCs
numSlicesSubUC = ceil(sDiff.numSlices/numSubUCs);
numSlices = sDiff.numSlices;

% store total intensity in the EW (optional for convergence testing)
% storeTotalInt = false

% coefficients
numUCs = coefs(1); % can be fractional
theta_x = coefs(2); %  hor tilt
theta_y = coefs(3); % vert tilt

sDiff.theta_x = theta_x;
sDiff.theta_y = theta_y;

% Compute potential at prescribed sample tilt
if sDiff.correctPotsForTilt
    sDiff = computePotential(sDiff);
    expPot = sDiff.expPot;
end

% Recompute propagator with tilts
% For uneven slices, sliceThickness is a vector with numSlices
% elements containing the thicknesses of each slice
t = reshape(sDiff.sliceThickness,...
    [1 1 length(sDiff.sliceThickness)]);
sDiff.prop = exp( ...
    sDiff.propConst ...
    + (2i*pi*t.* ...
    (tan(theta_x)*sDiff.qxa + tan(theta_y)*sDiff.qya)) ...
    );
if sDiff.flagUseAntiAliasing
    sDiff.prop = sDiff.qMask .* sDiff.prop;
end

% Determine if slices are empty or not (could do in setup phase)
isEmptySlice = false(sDiff.numSlices,1);
for iSlice = 1:numSlices
    isEmptySlice(iSlice) = ~any(reshape(sDiff.pot(:,:,iSlice),1,[]));
end

% Establish GPU arrays for computations
if nargin < 3
    expPot = sDiff.expPot;
    expPot = gpuArray(single(expPot)); % potentials too large to store on
    % GPU
end
prop = gpuArray(single(sDiff.prop));
EW = gpuArray(single(sDiff.EWInit));
%sigma = gpuArray(single(sDiff.sigma));

if flagStoreEWs
    EWstore = gpuArray(single(...
        zeros([sDiff.storeLenX,sDiff.storeLenY,round(numUCs*numSubUCs)])));
else
    EWstore = [];
end
% if storeTotalInt
%     IntStore = gpuArray(single(...
%         zeros(1,round(numUCs*numSubUCs))));
% else
%     IntStore = [];
% end

% measure time associated with various types of computation
% timeDiff = 0;
% timeStore = 0;

for a0 = 1:ceil(numUCs)
    if a0 == ceil(numUCs) % last 
        numSlices = round((numUCs-(a0-1))*sDiff.numSlices);
    end

    if computePot
        sDiff = computePotential(sDiff);
        expPot = gpuArray(single(exp(1i*sigma*sDiff.pot)));
    end
    
    for a1 = 1:numSlices

%         % Retrieve slice potential
%         if ndims(sDiff.pot) > 3
%             sliceExpPot = expPot(:,:,a1,randi(sDiff.numFP));
%         else
%             sliceExpPot = expPot(:,:,a1);
%         end
%         
%         % Retrieve slice propagator (needed for uneven slices)
%         if size(sDiff.prop,3) > 1
%             sliceProp = prop(:,:,a1);
%         else
%             sliceProp = prop;
%         end

%         tic
        if ~isEmptySlice(iSlice) % if potential slice is nonzero, interact and propagate
            EW(:) = prop(:,:,a1) .* fft2(ifft2(EW) .* ...
                    expPot(:,:,a1,randi(sDiff.numFP)));
        else % if potential slice is empty, just propagate
            EW(:) = sliceProp .* EW;
        end
%         tDiff = toc;
%         timeDiff = timeDiff + tDiff;
%      
%         tic
        % Store output
        if flagStoreEWs && mod(a1,numSlicesSubUC) == 0
            indStore = (a0-1)*numSubUCs+(a1/numSlicesSubUC);
            
            EWstore(:,:,indStore) = reshape(EW(sDiff.storeMask),...
                [sDiff.storeLenX,sDiff.storeLenY]);
%             if storeTotalInt
%                 IntStore(indStore) = sum(abs(EW(:)).^2); %% Istore (New)
%             end
        end
%         tStore = toc;
%         timeStore = timeStore + tStore;

    end
     
end

% disp(['Diffracting time (seconds): ' num2str(timeDiff)])
% disp(['Storing time (seconds): ' num2str(timeStore)])

end