%% Bloch Waves diffraction calculation

% Input variables
E0 = 750e3; %eV

% Relevant physical constants
lambElec = computeElectronWavelength(E0); % Electron wavelength (Angstroms)
m = 9.11e-31; % Electron mass (kg)
c = 3e8; % Speed of light (m/s)
e = 1.602e-19; % Charge of electron (C)
sigma = (2*pi/lambElec/E0) ...
    *(m*c^2+e*E0)/(2*m*c^2+e*E0); % Interaction parameter (rad / (V*A))
a0 = 0.5292; % Bohr radius in Angstroms

%% 1. Generate crystal structure and reciprocal lattice

% Generate crystal structure: unit cell and atomic coordinates
options.cellMult = 1;
[atoms,cellDim,lattice,uvwInit] = wyckoffGold(options);
Z = atoms(:,4); % Atomic number
nAtoms = size(lattice,1); % Number of atoms in unit cell
Vcryst = prod(cellDim); % unit cell volume (Angstroms^3)

% Compute reciprocal lattice vectors by real-space lattice inversion
Gvec = inv(uvwInit'); % Matrix of reciprocal lattice (row) vectors (Angstroms^-1)

% Generate mesh of hkl
hRange = [-26 26];
kRange = [-26 26];
lRange = [-2 2];
hLen = diff(hRange)+1;
kLen = diff(kRange)+1;
lLen = diff(lRange)+1;  
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

%% 2. Compute library of 3D specimen potential Fourier components U_G

% Compute structure factors

F_hkl = zeros(nOrders,1);
gamma = 1+(e*E0/m/c^2); % Electron Lorentz factor

scatApprox = 'Born'; % 'Born','Moliere'
% NOTE: Moliere approx gives complex scattering factors which makes
% U_G complex and in turn the matrix eigenvalues complex later on!
% Not sure if physical...
% NOTE: Check validity of the gamma correction to the
% scattering factors! Increases scattering factors for 
% relativistic beams... maybe not needed with sigma included?
correctScat = false;

for ii = 1:nOrders
    for xx = 1:nAtoms
        switch scatApprox
            case 'Born'
                if correctScat
                    f = gamma*electronScatteringFactor(Z(xx),Gmag(ii));
                else
                    f = electronScatteringFactor(Z(xx),Gmag(ii));
                end
            case 'Moliere'
                f = electronScatteringFactorMoliere(Z(xx),Gmag(ii),E0);
        end
        F_hkl(ii) = F_hkl(ii) ...
            + f *exp(-2i*pi*(hkl(ii,1)*lattice(xx,1) + hkl(ii,2)*lattice(xx,2)...
            + hkl(ii,3)*lattice(xx,3)));       
    end
end

U_G = (sigma/(pi*lambElec))*(47.86/Vcryst).*F_hkl;

% Debye-Waller factor
urms = 0.0894; % 1D rms displacement perpendicular to Bragg planes (Angstroms)
B = 8*pi^2*urms^2; % Debye-Waller B Factor
U_G = U_G.*exp(-B*Gmag.^2./4);

% Force U_G to be real (correct for centrosymmetric crystals)
U_G = real(U_G); 

%% 3. Compute excitation errors

% Sample tilt
thetaX = 0.1; % rad, x component of tilt
thetaY = 0; % rad, y component of tilt

s_G = computeExcitationError(thetaX,thetaY,Ghkl,lambElec); 
% Inverse Angstroms

%% 4. Select beams to include

U_thresh = 5e-3;
s_thresh = 5e-2;
isSel = abs(s_G) < s_thresh & abs(U_G) > U_thresh;
hklSel = hkl(isSel,:);
GhklSel = Ghkl(isSel,:);

%% 5. Build and solve matrix equation

N = sum(isSel);
k0 = 1./lambElec;
k0z = k0*(-cos(thetaX)*cos(thetaY));

hklDiff = zeros(N,N,3);
for ii = 1:3
    hklDiff(:,:,ii) = repmat(hklSel(:,ii),[1 N])...
        -repmat(hklSel(:,ii)',[N 1]);
end

indDiff = sub2ind([hLen,kLen,lLen],...
    hklDiff(:,:,1)-(hRange(1)-1),...
    hklDiff(:,:,2)-(kRange(1)-1),...
    hklDiff(:,:,3)-(lRange(1)-1));

A = diag(2*k0.*s_G(isSel))+U_G(indDiff).*(ones(N,N)-diag(ones(N,1)));

[Cvecs,eigvals] = eig(A);
Cvecsinv = conj(Cvecs');
% Cvecsinv = inv(Cvecs);

%% Initial condition
psiInit_G = zeros(N,1);
indZero = find(hklSel(:,1) == 0 ...
    & hklSel(:,2) == 0 ...
    & hklSel(:,3) == 0);
psiInit_G(indZero) = 1;

gam = eigvals/(2*k0z); % Inverse Angstroms

psi_G = @(z) Cvecs*(eye(N).*exp(2i*pi*gam*z))*Cvecsinv*psiInit_G;

%% 6. Generate thickness dependent scattered beams

nUC = 100;
dz = cellDim(3);
zTest = (0:nUC)*dz;
nZ = numel(zTest);

psi_G_array = zeros(nZ,N);

for iZ = 1:nZ
    psi_G_array(iZ,:) = psi_G(zTest(iZ));
end

figure;
plot(zTest./10,abs(psi_G_array).^2,'-')
xlabel('Distance (nm)')
ylabel('Beam amplitude')

%% 7. Generate thickness dependent diffraction pattern stack

Nx = 32;
Ny = Nx;
qx = makeFourierCoords(Nx,cellDim(1)/Nx);
qy = makeFourierCoords(Ny,cellDim(2)/Ny);
dqx = qx(2)-qx(1);
dqy = qy(2)-qy(1);
[qya,qxa] = meshgrid(qy,qx);
% q2 = sDiff.qxa.^2 + sDiff.qya.^2;

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

tArray = 0.1*zTest; % nm
StackViewerDiff(fftshift(fftshift(IDiff,1),2),tArray)
