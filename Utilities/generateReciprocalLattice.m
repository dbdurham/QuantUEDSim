function [hkl,Ghkl,Gmag,dhkl,Gvec] = generateReciprocalLattice(...
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

[Ghkl,Gmag,dhkl] = computeScatteringVectors(hkl,Gvec);

end

