function [Ilib,epsLib] = computeTiltAveragedDiffraction(...
    sigmaThetaSamp,nUC,nIter,diffMethod,sDiff,varargin)
%COMPUTETILTAVERAGEDDIFFRACTION Compute tilt-averaged diffraction pattern
%library using one of the available methods
%   sigmaThetaSamp = RMS tilt spreads to sample (rad)
%   nUC = Number of unit cells to sample
%   nIter = Compute up to this many iterations
%   diffMethod = 'Kinematical', 'Bloch Waves', or 'Multislice'
%   sDiff = struct containing inputs for the diffraction computation
%   function
%
%   Optional inputs:
%   Ilib = Previously started library of diffraction patterns
%   iStart = Iteration from which to continue the integration

useGPU = false;

sigmaThetaMax = max(sigmaThetaSamp);
nTheta = numel(sigmaThetaSamp);
symmDPs = true;

if nargin > 5
    % Extend computation of a previous library
    Ilib = varargin{1};
    iStart = varargin{2};
else
    % Initialize libraries
    Ilib = zeros([sDiff.storeSize,...
        nUC,...
        nTheta,...
        nIter]);
    iStart = 1;
end

if nargin > 7
    % extend computation of a subset of the tilt range
    tiltSubFactor = varargin{3};
else
    tiltSubFactor = 1;
end
indsThetaUpdate = 1:nTheta/tiltSubFactor; % all by default, subset if desired
indsThetaKeep = nTheta/tiltSubFactor+1:nTheta;

switch diffMethod
    case 'Kinematical'
        funcDiff = @(theta1,theta2) calcDiffKin(theta1,theta2,nUC,sDiff);        
    case 'Bloch Waves'
        funcDiff = @(theta1,theta2) calcDiffBW(theta1,theta2,nUC,sDiff);    
    case 'Multislice'
        funcDiff = @(theta1,theta2) calcDiffMS(theta1,theta2,nUC,sDiff,useGPU);
end

func = @(theta1,theta2) computeWeightedDiffStack(...
    theta1,theta2,sigmaThetaSamp(indsThetaUpdate),symmDPs,funcDiff);

iEnd = nIter;

% Only compute integral over one quadrant if applying four-fold symmetry
if symmDPs
    lowerTheta = 0;
else
    lowerTheta = -3*sigmaThetaMax/tiltSubFactor;
end
upperTheta = 3*sigmaThetaMax/tiltSubFactor;

iterSub = log2(tiltSubFactor);

for iIter = iStart:iEnd
    tic
    if iIter == 1
        Ilib(:,:,:,indsThetaUpdate,1) = extendTrapz2D(func,...
            lowerTheta,upperTheta,...
            lowerTheta,upperTheta,...
            Ilib(:,:,:,indsThetaUpdate,1),...
            1);
    else
        Ilib(:,:,:,indsThetaUpdate,iIter) = extendTrapz2D(func,...
            lowerTheta,upperTheta,...
            lowerTheta,upperTheta,...
            Ilib(:,:,:,indsThetaUpdate,iIter-1),...
            iIter-iterSub);
        Ilib(:,:,:,indsThetaKeep,iIter) = ...
            Ilib(:,:,:,indsThetaKeep,iIter-1);
    end
    
    
    if symmDPs
        % Apply vertical mirror
        Istack2 = (Ilib(:,:,:,:,iIter)...
            + Ilib([1 end:-1:2],:,:,:,iIter))./2;
        % Apply horizontal mirror
        Istack2 = (Istack2...
            + Istack2(:,[1 end:-1:2],:,:))./2;
        % Apply 90 deg rotation
        rotStack = rot90(fftshift(fftshift(Istack2,1),2));
        if mod(sDiff.storeSize(1),2) == 0
            rotStack = circshift(rotStack,1,1);
        end
        Ilib(:,:,:,:,iIter) = ifftshift(ifftshift(...
            (fftshift(fftshift(Istack2,1),2) ...
            + rotStack)./2 ...
            ,1),2);
    end
    clear('Istack2')
        
    disp(['Computed iteration #' num2str(iIter)])
    toc
    
    if iIter > 1
        % Fractional change as function of thickness and tilt spread
        epsLib = squeeze(sum(sum(abs(...
            Ilib(:,:,:,indsThetaUpdate(end),iIter) - Ilib(:,:,:,indsThetaUpdate(end),iIter-1)),1),2)...
            ./ sum(sum(Ilib(:,:,:,indsThetaUpdate(end),iIter-1),1),2));
        disp(['Fractional change in diff stack at max spread: ' num2str(mean(epsLib))])
    end
    
end

end

