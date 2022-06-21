function [Iarray,psi_G_array,hklSel,GhklSel] = calcIntsBW(theta1,theta2,...
    nUC,sDiff)
%CALCDIFFBW Calculate diffracted intensities using Bloch Wave method
%   theta1 = x component of sample tilt (rad)
%   theta2 = y component (rad)
%   nUC = number of unit cells to sample
%   sDiff = setup struct
%   NOTE: ASSUMES CENTROSYMMETRIC CRYSTAL

% Compute excitation errors (inv Angstroms)
s_G = computeExcitationError(theta1,theta2,sDiff.Ghkl,sDiff.lambElec); 

% Select beams within specified reciprocal space volume around 
% the Ewald sphere
Gxy = sqrt(sDiff.Ghkl(:,1).^2 + sDiff.Ghkl(:,2).^2);
% Gz = abs(Ghkl(:,3));
isSel = Gxy < sDiff.GxyThresh ...
    & ~(sDiff.U_G==0) ...
    & abs(s_G) < sDiff.sThresh;
hklSel = sDiff.hkl(isSel,:);
GhklSel = sDiff.Ghkl(isSel,:);


%% Build and solve matrix equation

N = sum(isSel);
k0 = 1./sDiff.lambElec;
k0z = k0*(-cos(theta1)*cos(theta2));

hklDiff = zeros(N,N,3);
for ii = 1:3
    hklDiff(:,:,ii) = repmat(hklSel(:,ii),[1 N])...
        -repmat(hklSel(:,ii)',[N 1]);
end

if any(any(hklDiff(:,:,1) > sDiff.hRange(2) ...
        | hklDiff(:,:,1) < sDiff.hRange(1)))
    disp('h out of range, need more h points')
end
if any(any(hklDiff(:,:,2) > sDiff.kRange(2) ...
        | hklDiff(:,:,2) < sDiff.kRange(1)))
    disp('k out of range, need more k points')
end
if any(any(hklDiff(:,:,3) > sDiff.lRange(2) ...
        | hklDiff(:,:,3) < sDiff.lRange(1)))
    disp('l out of range, need more l points')
end


indDiff = sub2ind([sDiff.hLen,sDiff.kLen,sDiff.lLen],...
    hklDiff(:,:,1)-(sDiff.hRange(1)-1),...
    hklDiff(:,:,2)-(sDiff.kRange(1)-1),...
    hklDiff(:,:,3)-(sDiff.lRange(1)-1));

A = diag(2*k0.*s_G(isSel))+sDiff.U_G(indDiff).*(ones(N,N)-diag(ones(N,1)));

[Cvecs,eigvals] = eig(A);
Cvecsinv = conj(Cvecs');


% Cvecsinv = inv(Cvecs);

%% Compute exit wave components
% Initial condition (plane wave)
psiInit_G = zeros(N,1);
psiInit_G(hklSel(:,1)==0 & hklSel(:,2)==0 & hklSel(:,3)==0) = 1;

gamVec = diag(eigvals)/(2*k0z); % Inverse Angstroms
RHS = Cvecsinv*psiInit_G;
psi_G = @(z) Cvecs*(exp(2i*pi*gamVec*z).*RHS);

dz = sDiff.cellDim(3);
zTest = (1:nUC)*dz;
nZ = numel(zTest);

psi_G_array = zeros(N,nZ);

for iZ = 1:nZ
    psi_G_array(:,iZ) = psi_G(zTest(iZ));
end

Iarray = abs(psi_G_array).^2;

end

