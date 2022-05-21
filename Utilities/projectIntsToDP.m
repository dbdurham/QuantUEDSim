function [IDiff,qxa,qya] = projectIntsToDP(IArray,GhklSel,sDiff)
%PROJECTINTsTODP Calculate diffraction pattern from Bloch wave computed
%intensities
%   Iarray = Diffraction peak intensities vs thickness
%   Ghklsel = Reciprocal lattice vectors of the selected peaks
%   sDiff - struct containing sim parameters
%

N = size(IArray,1);
nUC = size(IArray,2);

% Unpack input variable structure
fieldNames = fieldnames(sDiff);
nFields = numel(fieldNames);
for iField = 1:nFields
    [~] = evalc([fieldNames{iField} ' = sDiff.' fieldNames{iField}]);
end

% Generate reciprocal-space grid on which to project DP
Nx = 32;
Ny = Nx;
qx = makeFourierCoords(Nx,cellDim(1)/Nx);
qy = makeFourierCoords(Ny,cellDim(2)/Ny);
[qya,qxa] = meshgrid(qy,qx);

% Planar approximation to find diffraction positions
GmagSel = sqrt(GhklSel(:,1).^2 + GhklSel(:,2).^2 + GhklSel(:,3).^2);
GthetaSel = atan2(GhklSel(:,2),GhklSel(:,1));
GxProj = GmagSel.*cos(GthetaSel);
GyProj = GmagSel.*sin(GthetaSel);
% identify indices to fill in patterns
indxProj = zeros(N,1);
indyProj = zeros(N,1);
for iBeam = 1:N
    [~,indxProj(iBeam)] = min(abs(GxProj(iBeam)-qx));
    [~,indyProj(iBeam)] = min(abs(GyProj(iBeam)-qy));
end

% Accumulte diffracted intensities on the generated grid
IDiff = zeros(Ny,Nx,nUC);
for iUC = 1:nUC
    IDiff(:,:,iUC) = accumarray([indxProj,indyProj],...
        IArray(:,iUC)',[Ny,Nx]);
end

end

