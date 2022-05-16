function [Istack,weights,inSim] = computeWeightedIntsKin(theta1,theta2,...
    sigmaThetaSamp,sigmaThetaMax,t,Ghkl,SFMag,lamElec,Vcryst,alpha,B,symmDPs)
%COMPUTEWEIGHTEDINTSKIN
%   

thetasqr = theta1.^2 + theta2.^2;
weights = (1./(2*pi*sigmaThetaSamp.^2)).*exp(-thetasqr./(2*sigmaThetaSamp.^2)) ... % Gaussian weighting
    .* (abs(theta1).^2 + abs(theta2).^2 <= (3*sigmaThetaSamp).^2) ... % radial mask up to 3 sigma
    .* 0.5.^((abs(theta1) == 3*sigmaThetaSamp).*(sigmaThetaSamp<sigmaThetaMax)) ... % half weight along x edge (trapezoid rule)
    .* 0.5.^((abs(theta2) == 3*sigmaThetaSamp).*(sigmaThetaSamp<sigmaThetaMax)); % half weight along y edge (trapezoid rule)
if symmDPs 
   weights = weights .* 4 ... % correct for only using 1/4 of the region
       .* 0.5.^((abs(theta1) == 0).*(sigmaThetaSamp<sigmaThetaMax)) ... % half weight along x axis
       .* 0.5.^((abs(theta2) == 0).*(sigmaThetaSamp<sigmaThetaMax)); % half weight along y axis
end
nTheta = numel(weights);
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


