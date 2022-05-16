function [I,excitationError] = computeDiffKin(t,thetaX,thetaY,...
    Ghkl,SFMag,lamElec,Vcryst,alpha,B)
%COMPUTEDIFFKIN Compute kinematical diffraction pattern
%   t = thickness (Angstroms), can be a row vector
%   thetaX, thetaY = sample tilt angles in rad
%   Ghkl = N x 3 array of g vectors (Angstroms^-1)
%   braggAngle = N x 1 array of bragg angles
%   SFMag = N x 1 array of structure factor magnitudes
%   lamElec = electron de Broglie wavelength (Angstroms)
%   Vcryst = unit cell volume (Angstroms^3)
%   alpha = absorption term
%   B = Debye-Waller factor

nOrders = size(Ghkl,1);
Gmag = sqrt(Ghkl(:,1).^2 + Ghkl(:,2).^2 + Ghkl(:,3).^2);
dhkl = 1./Gmag;
braggAngle = asin(lamElec./(2*dhkl)); % Bragg angle (rad)

% Incident beam vector
% NOTE: ASSUMES 001 ZONE AXIS, NEED TO GENERALIZE FOR ARBITRARY ZONE
% zoneAxis = [0; 0; 1]; % Zone axis given in miller indices
% zoneCart = sum(zoneAxis.*uvwInit,1); 
% zoneCart = zoneCart./norm(zoneCart); % Zone axis direction in cartesian coords
kiMag = 1./lamElec; % Incident electron wave vector magnitude (Angstroms ^ -1)
ki = kiMag.*[cos(thetaY)*sin(thetaX) ...
    cos(thetaX)*sin(thetaY) ...
    -cos(thetaX)*cos(thetaY)]; % Wave vector
betaAngle = acos(cos(thetaX)*cos(thetaY)); % Angle between incident beam and surface normal (rad)
% (0 at normal incidence)

extDist = pi.*Vcryst.*cos(braggAngle)... % Extinction distance (Angstroms)
    ./(lamElec.*SFMag); 
excitationError = zeros(nOrders,1); % Excitation error (Angstroms^-1)

for ii = 1:nOrders
    excitationError(ii)  = -dot(Ghkl(ii,:),2*ki+Ghkl(ii,:))...
        ./(2*norm(ki+Ghkl(ii,:))*cos(betaAngle));
end

I = (sin(pi*excitationError*t)./(excitationError.*extDist)).^2 ...
    .*exp(-alpha*t).*exp(-B*Gmag.^2./2); % Diffracted intensities vs thickness
% Handle zero-order beam (excitation error = 0)
for ii = 1:nOrders
    if excitationError(ii) == 0
        I(ii,:) = (pi*t./extDist(ii)).^2 ...
            .*exp(-alpha*t); % Diffracted intensities vs thickness
    end
end

end

