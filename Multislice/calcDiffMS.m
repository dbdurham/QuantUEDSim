function IDiff = calcDiffMS(theta1,theta2,numUC,sDiff)
%CALCDIFFMS Wrapper function for multislice diffraction calculation
%   

if sDiff.useGPU
    [~,EWStore,~] = calcDiffMSGPU(sDiff,...
        [numUC,theta1,theta2],sDiff.expPot);
    IDiff = double(abs(gather(EWStore)).^2);
else
    [~,EWStore,~] = calcDiffMSCPU(sDiff,...
        [numUC,theta1,theta2]);
    IDiff = abs(EWStore).^2;
end

end

