function [IDiff,qxa,qya] = calcDiffBW(psi_G_array,isSel,sDiff)
%CALCDIFFBW Calculate diffraction pattern from Bloch wave computed
%intensities
%   psi_G_Array = Bloch wave computed intensities
%   sDiff
%   NOTE: ASSUMES CENTROSYMMETRIC CRYSTAL

N = sum(isSel);
nUC = size(psi_G_array,1);

% Unpack input variable structure
fieldNames = fieldnames(sDiff);
nFields = numel(fieldNames);
for iField = 1:nFields
    [~] = evalc([fieldNames{iField} ' = sDiff.' fieldNames{iField}]);
end

Nx = 32;
Ny = Nx;
qx = makeFourierCoords(Nx,cellDim(1)/Nx);
qy = makeFourierCoords(Ny,cellDim(2)/Ny);
dqx = qx(2)-qx(1);
dqy = qy(2)-qy(1);
[qya,qxa] = meshgrid(qy,qx);

% Planar approximation to find diffraction positions
Gtheta = atan2(Ghkl(:,2),Ghkl(:,1));
GxProj = Gmag(isSel).*cos(Gtheta(isSel));
GyProj = Gmag(isSel).*sin(Gtheta(isSel));
% identify indices to fill in patterns
indxProj = zeros(N,1);
indyProj = zeros(N,1);
for iBeam = 1:N
    [~,indxProj(iBeam)] = min(abs(GxProj(iBeam)-qx));
    [~,indyProj(iBeam)] = min(abs(GyProj(iBeam)-qy));
end

IDiff = zeros(Ny,Nx,nUC);
for iUC = 1:nUC
    IDiff(:,:,iUC) = accumarray([indxProj,indyProj],...
        abs(psi_G_array(iUC,:)').^2,[Ny,Nx]);
end

end

