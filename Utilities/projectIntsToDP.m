function IDiff = projectIntsToDP(IArray,GhklSel,qxa,qya)
%PROJECTINTsTODP Calculate diffraction pattern from Bloch wave computed
%intensities
%   Iarray = Diffraction peak intensities vs thickness
%   Ghklsel = Reciprocal lattice vectors of the selected peaks
%   qxa - reciprocal space coord matrix (x)
%   qya - likewise for y
%   NOTE: qx goes along rows (1), qy along columns (2)

N = size(IArray,1);
nUC = size(IArray,2);

% Planar approximation to find diffraction positions
GmagSel = sqrt(GhklSel(:,1).^2 + GhklSel(:,2).^2 + GhklSel(:,3).^2);
GthetaSel = atan2(GhklSel(:,2),GhklSel(:,1));
GxProj = GmagSel.*cos(GthetaSel);
GyProj = GmagSel.*sin(GthetaSel);
% identify indices to fill in patterns
indxProj = zeros(N,1);
indyProj = zeros(N,1);
qx = qxa(:,1);
qy = qya(1,:);
for iBeam = 1:N
    [~,indxProj(iBeam)] = min(abs(GxProj(iBeam)-qx));
    [~,indyProj(iBeam)] = min(abs(GyProj(iBeam)-qy));
end

% Accumulte diffracted intensities on the generated grid
Nx = numel(qx);
Ny = numel(qy);
IDiff = zeros(Nx,Ny,nUC);
for iUC = 1:nUC
    IDiff(:,:,iUC) = accumarray([indxProj,indyProj],...
        IArray(:,iUC)',[Nx,Ny]);
end

end

