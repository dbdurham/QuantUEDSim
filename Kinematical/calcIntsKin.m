function I = calcIntsKin(theta1,theta2,nUC,sDiff)
%CALCINTSKIN Calculate diffracted intensities using kinematical model
%   theta1 = x component of sample tilt (rad)
%   theta2 = y component (rad)
%   sDiff = setup struct

t = (1:nUC)*sDiff.cellDim(3);

% Compute excitation errors (inv Angstroms)

s_G = computeExcitationError(theta1,theta2,sDiff.Ghkl,sDiff.lambElec); 


% Calculate extinction distances (Angstroms)
braggAngle = asin(sDiff.lambElec./(2*sDiff.dhkl)); % Bragg angle (rad)
extDist = pi.*sDiff.Vcryst.*cos(braggAngle)...
    ./(sDiff.lambElec.*abs(sDiff.Fhkl));

% Calculate diffracted intensities vs thickness
I = (sin(pi*s_G*t)./(s_G.*extDist)).^2;
% Handle zero-order beam (excitation error = 0)
I(s_G==0,:) = (pi*t./extDist(s_G==0)).^2;
% Apply DWF
[~,DWInt] = computeDWF(sDiff.uRMS,1,sDiff.Gmag);
I = I.*DWInt;

end

