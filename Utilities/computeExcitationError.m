function [excitationError] = computeExcitationError(...
    thetaX,thetaY,Ghkl,lambElec)
%COMPUTEEXCITATIONERROR Computes the excitation error in inverse 
% Angstroms for a reciprocal lattice illuminated by radiation 
% with a given sample orientation and wavelength
%   thetaX, thetaY = beam-crystal angle
%   Ghkl = reciprocal lattice vector list (N x 3 array)
%   lambElec = wavelength of incident radiation

kiMag = 1./lambElec; % Incident electron wave vector magnitude (Angstroms ^ -1)
ki = kiMag.*[cos(thetaY)*sin(thetaX) ...
    cos(thetaX)*sin(thetaY) ...
    -cos(thetaX)*cos(thetaY)]; % Wave vector
betaAngle = acos(cos(thetaX)*cos(thetaY)); % Angle between incident beam and surface normal (rad)
% (0 at normal incidence)

nOrders = size(Ghkl,1);
excitationError = zeros(nOrders,1); % Excitation error (Angstroms^-1)
for ii = 1:nOrders
    excitationError(ii)  = -dot(Ghkl(ii,:),2*ki+Ghkl(ii,:))...
        ./(2*norm(ki+Ghkl(ii,:))*cos(betaAngle));
end

end

