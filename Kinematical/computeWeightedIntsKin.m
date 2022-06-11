function [Istack,weights,inSim] = computeWeightedIntsKin(theta1,theta2,...
    sigmaThetaSamp,sigmaThetaMax,t,Ghkl,SFMag,lamElec,Vcryst,alpha,B,symmDPs)
%COMPUTEWEIGHTEDINTSKIN
%   

weights = computeGaussianTiltSpreadWeights(theta1,theta2,...
    sigmaThetaSamp,symmDPs);

nTheta = numel(sigmaThetaSamp);
inSim = (abs(theta1).^2 + abs(theta2).^2 <= (3*sigmaThetaMax)^2);

if inSim
    I = computeDiffKin(...
        t,... % thicknesses
        theta1,theta2,... % sample tilt angles
        Ghkl,SFMag,... % g vectors and structure factors of reflections
        lamElec,... % electron wavelength
        Vcryst,... % unit cell volume
        alpha,... % absorption length
        B); % Debye-Waller Factor
    Istack = zeros([size(I),nTheta]);
    for iTheta = 1:nTheta
        Istack(:,:,iTheta) = I .* weights(iTheta);
    end
else
    Istack = 0;
end

end


