function weightedDiffStack = computeWeightedDiffStack(theta1,theta2,...
    sigmaThetaSamp,symmDPs,funcDiff)
%COMPUTEWEIGHTEDDIFFSTACK Summary of this function goes here
%   Detailed explanation goes here
nTheta = numel(sigmaThetaSamp);
diffMat = repmat(funcDiff(theta1,theta2),[1 1 1 nTheta]);
diffSize = size(diffMat);
weightMat = repmat(permute(computeGaussianTiltSpreadWeights(theta1,theta2,...
            sigmaThetaSamp,symmDPs),[3 4 1 2]),[diffSize(1:2) 1]);
weightedDiffStack = weightMat.*diffMat;

end

