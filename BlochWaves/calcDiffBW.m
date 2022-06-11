function IDiffStore = calcDiffBW(theta1,theta2,nUC,sDiff)
%CALCDIFFBW Calculate diffraction pattern stack using Bloch Wave method
%   Detailed explanation goes here

[IArray,~,~,GhklSel] = calcIntsBW(theta1,theta2,nUC,sDiff);
IDiff = projectIntsToDP(IArray,GhklSel,sDiff.qxa,sDiff.qya);
% Downsample for storage
IDiffStore = IDiff(1:sDiff.downSampFac:end,1:sDiff.downSampFac:end,:);

end

