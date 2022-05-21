function IArray = extractsIntsFromDP(IDiff,qxa,qya,...
    hklExt,sDiff)
%EXTRACTSINTSFROMDP Retrieve diffracted intensities from a simulated
%diffraction pattern stack at specified orders
%   IDiff - MxNxP stack of simulated diffraction patterns
%   qxa - MxN reciprocal space coord array (x)
%   qxy - likewise for y
%   hklSel - Qx3 selected miller indices at which to extract intensities
%   sDiff - struct containing sim parameters
%   NOTE: code assumes that specified diffraction orders appear in the
%   pattern

N = size(hklExt,1);
nUC = size(IDiff,3);

% Unpack input variable structure
fieldNames = fieldnames(sDiff);
nFields = numel(fieldNames);
for iField = 1:nFields
    [~] = evalc([fieldNames{iField} ' = sDiff.' fieldNames{iField}]);
end

qx = qxa(1,:);
qy = qya(:,1);

% Compute scattering vectors
[GhklSel,GmagSel] = computeScatteringVectors(hkl,Gvec);

% Planar approximation to find diffraction positions
GthetaSel = atan2(GhklSel(:,2),GhklSel(:,1));
GxProj = GmagSel.*cos(GthetaSel);
GyProj = GmagSel.*sin(GthetaSel);
% identify indices to extract from patterns
indxProj = zeros(N,1);
indyProj = zeros(N,1);
for iBeam = 1:N
    [~,indxProj(iBeam)] = min(abs(GxProj(iBeam)-qx));
    [~,indyProj(iBeam)] = min(abs(GyProj(iBeam)-qy));
end

% Extract intensities
IArray = zeros(N,nUC);
for iBeam = 1:N
    IArray(iBeam,:) = IDiff(indyProj(iBeam),indxProj(iBeam),:);
end

end

