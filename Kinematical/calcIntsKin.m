function I = calcIntsKin(theta1,theta2,nUC,sDiff)
%CALCINTSKIN Calculate diffracted intensities using kinematical model
%   theta1 = x component of sample tilt (rad)
%   theta2 = y component (rad)
%   sDiff = setup struct

% Unpack input variable structure
fieldNames = fieldnames(sDiff);
nFields = numel(fieldNames);
for iField = 1:nFields
    [~] = evalc([fieldNames{iField} ' = sDiff.' fieldNames{iField}]);
end

t = (1:nUC)*cellDim(3);

% Compute excitation errors (inv Angstroms)
s_G = computeExcitationError(theta1,theta2,Ghkl,lambElec); 

% Calculate extinction distances (Angstroms)
braggAngle = asin(lambElec./(2*dhkl)); % Bragg angle (rad)
extDist = pi.*Vcryst.*cos(braggAngle)./(lambElec.*abs(Fhkl));

% Calculate diffracted intensities vs thickness
I = (sin(pi*s_G*t)./(s_G.*extDist)).^2;
% Handle zero-order beam (excitation error = 0)
I(s_G==0,:) = (pi*t./extDist(s_G==0)).^2;
% Apply DWF
[~,DWInt] = computeDWF(uRMS,1,Gmag);
I = I.*DWInt;

end

