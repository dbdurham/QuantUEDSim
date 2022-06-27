function IDiff = calcDiffMS(theta1,theta2,numUC,sDiff,useGPU,varargin)
%CALCDIFFMS Wrapper function for multislice diffraction calculation
%   

if nargin > 5
    expPot = varargin{1};
end

if useGPU
    [~,EWStore,~] = calcDiffGPU(sDiff,...
        [numUC,0,theta1,theta2],...
        expPot);
    IDiff = double(abs(gather(EWStore)).^2);
else
    [~,EWStore,~] = calcDiff(sDiff,...
        [numUC,0,theta1,theta2]);
    IDiff = abs(EWStore).^2;
end

end

