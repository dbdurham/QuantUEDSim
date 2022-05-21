function [Ghkl,Gmag,dhkl] = computeScatteringVectors(hkl,Gvec)
%COMPUTESCATTERINGVECTORS Compute scattering vectors for a specified set of
%miller indices and reciprocal lattice vectors
%   hkl - N x 3 array of miller indices
%   Gvec - Matrix of reciprocal lattice (row) vectors
%   

nOrders = size(hkl,1);
Ghkl = zeros(nOrders,3);
Gmag = zeros(nOrders,1);
for iOrder = 1:nOrders
    Ghkl(iOrder,:) = hkl(iOrder,1)*Gvec(1,:) ...
        + hkl(iOrder,2)*Gvec(2,:) ...
        + hkl(iOrder,3)*Gvec(3,:) ;
    Gmag(iOrder) = norm(Ghkl(iOrder,:));
end
dhkl = 1./Gmag; % Interplanar spacing

end

