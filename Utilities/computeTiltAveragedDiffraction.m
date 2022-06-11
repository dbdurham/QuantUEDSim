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

sigmaThetaMax = max(sigmaThetaSamp);
symmDPs = true;

if nargin > 5
    % Extend computation of a previous library
    Ilib = varargin{1};
    iStart = varargin{2};
else
    % Initialize libraries
    nTheta = numel(sigmaThetaSamp);
    Ilib = zeros([sDiff.storeSize,...
        nUC,...
        nTheta,...
        nIter]);
    iStart = 1;
end

switch diffMethod
    case 'Kinematical'
        funcDiff = @(theta1,theta2) calcDiffKin(theta1,theta2,nUC,sDiff);        
    case 'Bloch Waves'
        funcDiff = @(theta1,theta2) calcDiffBW(theta1,theta2,nUC,sDiff);    
    case 'Multislice'
        funcDiff = @(theta1,theta2) calcDiffMS(sDiff,[nUC,0,theta1,theta2]);
end

func = @(theta1,theta2) computeWeightedDiffStack(...
    theta1,theta2,sigmaThetaSamp,symmDPs,funcDiff);

iEnd = nIter;

% Only compute integral over one quadrant if applying four-fold symmetry
if symmDPs
    lowerTheta = 0;
else
    lowerTheta = -3*sigmaThetaMax;
end

for iIter = iStart:iEnd
    tic
    if iIter == 1
        Ilib(:,:,:,:,1) = extendTrapz2D(func,lowerTheta,3*sigmaThetaMax,...
            lowerTheta,3*sigmaThetaMax,Ilib(:,:,:,:,1),iIter);
    else
        Ilib(:,:,:,:,iIter) = extendTrapz2D(func,lowerTheta,3*sigmaThetaMax,...
            lowerTheta,3*sigmaThetaMax,Ilib(:,:,:,:,iIter-1),iIter);
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
        epsLib = squeeze(sum(abs(...
            Ilib(:,:,:,:,iIter) - Ilib(:,:,:,:,iIter-1)),[1 2])...
            ./ sum(Ilib(:,:,:,:,iIter-1),[1 2]));
        disp(['Fractional change in diff stack at max spread: ' num2str(mean(epsLib(:,end)))])
    end
    
end

end

