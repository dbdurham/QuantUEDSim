function IDiffStore = calcDiffKin(theta1,theta2,nUC,sDiff)
%CALCDIFFBW Calculate diffraction pattern stack using Kinematical method
%   Detailed explanation goes here

IArray = calcIntsKin(theta1,theta2,nUC,sDiff);
IDiff = projectIntsToDP(IArray,sDiff.Ghkl,sDiff.qxa,sDiff.qya);
% Downsample for storage
IDiffStore = IDiff(1:sDiff.downSampFac:end,1:sDiff.downSampFac:end,:);

end

