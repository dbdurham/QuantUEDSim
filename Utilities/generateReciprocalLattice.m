function [hkl,Ghkl,Gmag,dhkl] = generateReciprocalLattice(...
    uvwInit,hRange,kRange,lRange)
%GENERATERECIPROCALLATTICE Generate a list of reciprocal lattice vectors
%over a prescribed h, k, and l range
%   uvwInit: A 3 x 3 matrix with the real-space lattice vectors along the
%   columns (?)
%   hRange: 2 x 1 vector [lower_limit upper_limit] for h miller index
%   kRange: likewise for k
%   lRange: likewise for l

% Compute reciprocal lattice vectors by real-space lattice inversion
Gvec = inv(uvwInit'); % Matrix of reciprocal lattice (row) vectors (Angstroms^-1)

[h,k,l] = ndgrid(hRange(1):hRange(2),kRange(1):kRange(2),lRange(1):lRange(2));
hkl = [h(:) k(:) l(:)];
nOrders = size(hkl,1);

% Compute scattering vectors and magnitudes (Angstroms^-1)
Ghkl = zeros(nOrders,3);
Gmag = zeros(nOrders,1);
for iOrder = 1:nOrders
    Ghkl(iOrder,:) = hkl(iOrder,1)*Gvec(1,:) ...
        + hkl(iOrder,2)*Gvec(2,:) ...
        + hkl(iOrder,3)*Gvec(3,:) ;
    Gmag(iOrder) = norm(Ghkl(iOrder,:));
end
dhkl = 1./Gmag; % Interplanar spacing (Angstroms)

end

