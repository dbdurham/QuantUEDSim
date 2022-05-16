function Fhkl = computeStructureFactors(lattice,Z,hkl,Gmag,E0,...
    scatApprox,correctLorentz)
%COMPUTESTRUCTUREFACTORS Compute structure factor (in Angstroms) 
%   lattice - Nx3 array of relative coordinates of atoms in the unit cell
%   (from 0 to 1)
%   Z - Nx1 list of atomic numbers
%   hkl - Px3 array of h,k,l miller indices of desired reflections
%   Gmag - Px1 array of reciprocal lattice magnitudes
%   E0 - Kinetic energy of electrons (eV)
%   scatApprox - 'Born' or 'Moliere'
%   correctLorentz - Boolean option to correct Born approximation
%   atomic scattering factors by gamma, sometimes used for kinematical calcs of relativistic beams

load('fparams.mat','fparams')

nAtoms = size(lattice,1);
nOrders = size(hkl,1);
Fhkl = zeros(nOrders,1);

for ii = 1:nOrders
    for xx = 1:nAtoms
        switch scatApprox
            case 'Born'
                f = electronScatteringFactor(Z(xx),Gmag(ii),fparams);
                if correctLorentz
                    f = gamma*f;
                end
            case 'Moliere'
                f = electronScatteringFactorMoliere(Z(xx),Gmag(ii),E0,fparams);
        end
        Fhkl(ii) = Fhkl(ii) ...
            + f *exp(-2i*pi*(hkl(ii,1)*lattice(xx,1) + hkl(ii,2)*lattice(xx,2)...
            + hkl(ii,3)*lattice(xx,3)));    
    end
end

end

