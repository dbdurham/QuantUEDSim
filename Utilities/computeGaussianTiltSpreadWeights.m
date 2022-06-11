function weights = computeGaussianTiltSpreadWeights(theta1,theta2,...
    sigmaThetaSamp,symmDPs)
%COMPUTETILTSPREADWEIGHTS Compute the weights from the Gaussian tilt
%distribution at the sampled tilt angle(s)
%   theta1 = Tilt in x (rad)
%   theta2 = Tilt in y (rad)
%   sigmaThetaSamp = Tilt spreads being sampled
%   symmDPs = Boolean toggle for 4-fold symmetrization of patterns

sigmaThetaMax = max(sigmaThetaSamp);

thetasqr = theta1.^2 + theta2.^2;
weights = (1./(2*pi*sigmaThetaSamp.^2)).*exp(-thetasqr./(2*sigmaThetaSamp.^2)) ... % Gaussian weighting
    .* (abs(theta1).^2 + abs(theta2).^2 <= (3*sigmaThetaSamp).^2) ... % radial mask up to 3 sigma
    .* 0.5.^((abs(theta1) == 3*sigmaThetaSamp).*(sigmaThetaSamp<sigmaThetaMax)) ... % half weight along x edge (trapezoid rule)
    .* 0.5.^((abs(theta2) == 3*sigmaThetaSamp).*(sigmaThetaSamp<sigmaThetaMax)); % half weight along y edge (trapezoid rule)
if symmDPs 
   weights = weights .* 4; % correct for only using 1/4 of the region
       % .* 0.5.^((abs(theta1) == 0).*(sigmaThetaSamp<sigmaThetaMax)) ... % half weight along x axis
       % .* 0.5.^((abs(theta2) == 0).*(sigmaThetaSamp<sigmaThetaMax)); % half weight along y axis
end

end

