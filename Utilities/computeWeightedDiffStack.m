function weightedDiffStack = computeWeightedDiffStack(theta1,theta2,...
    sigmaThetaSamp,symmDPs,funcDiff)
%COMPUTEWEIGHTEDDIFFSTACK Summary of this function goes here
%   Detailed explanation goes here

diffMat = funcDiff(theta1,theta2);

weightMat = permute(computeGaussianTiltSpreadWeights(theta1,theta2,...
            sigmaThetaSamp,symmDPs),[3 4 1 2]);
weightedDiffStack = weightMat.*diffMat;

end

