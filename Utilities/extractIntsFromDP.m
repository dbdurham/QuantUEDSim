function IArray = extractIntsFromDP(IDiff,qxa,qya,GhklSel)
%EXTRACTSINTSFROMDP Retrieve diffracted intensities from a simulated
%diffraction pattern stack at specified orders
%   IDiff - MxNxP stack of simulated diffraction patterns
%   qxa - matrix of reciprocal space positions along x
%   qya - likewise for y
%   GhklSel - Q x 3 3D reciprocal space points to sample
%   NOTE: code assumes that specified diffraction orders appear in the
%   pattern

N = size(GhklSel,1);
nUC = size(IDiff,3);

% Planar approximation to find diffraction positions
GmagSel = sqrt(GhklSel(:,1).^2 + GhklSel(:,2).^2 + GhklSel(:,3).^2);
GthetaSel = atan2(GhklSel(:,2),GhklSel(:,1));
GxProj = GmagSel.*cos(GthetaSel);
GyProj = GmagSel.*sin(GthetaSel);
% identify indices to extract from patterns
indxProj = zeros(N,1);
indyProj = zeros(N,1);
qx = qxa(:,1);
qy = qya(1,:);
for iBeam = 1:N
    [~,indxProj(iBeam)] = min(abs(GxProj(iBeam)-qx));
    [~,indyProj(iBeam)] = min(abs(GyProj(iBeam)-qy));
end

% Extract intensities
IArray = zeros(N,nUC);
for iBeam = 1:N
    IArray(iBeam,:) = IDiff(indxProj(iBeam),indyProj(iBeam),:);
end

end

