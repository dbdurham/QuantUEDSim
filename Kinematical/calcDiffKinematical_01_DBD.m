clear all; clc; close all; 

%% Generate crystal structure and reciprocal lattice

% Generate crystal structure: unit cell and atomic coordinates
options.cellMult = 1;
[atoms,cellDim,lattice,uvwInit] = wyckoffGold(options);
Z = atoms(:,4); % Atomic number
nAtoms = size(lattice,1); % Number of atoms in unit cell

% Compute reciprocal lattice vectors by real-space lattice inversion
Gvec = inv(uvwInit'); % Matrix of reciprocal lattice (row) vectors (Angstroms^-1)

% Generate mesh of hkl
hRange = [-8 8];
kRange = [-8 8];
lRange = [0 0];
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

%% Calculate structure factors
% NOTE: scattering factors assumed independent of beam energy... an
% approximation often made for relativistic electrons

% Beam characteristics
E0 = 3000; %keV
lamElec = wavelengthElectrons(E0); % Electron de Broglie wavelength (Angstroms)
m = 9.109383*10^-31;
e = 1.602177*10^-19;
c =  299792458;
gamma = 1+(e*(E0*1e3)/m/c^2); % Electron Lorentz factor

F_hkl = zeros(nOrders,1); % Complex structure factor

scatApprox = 'Moliere'; % 'Born','Moliere'

for ii = 1:nOrders
    for xx = 1:nAtoms
        switch scatApprox
            case 'Born'
                f = gamma*electronScatteringFactor(Z(xx),Gmag(ii));
            case 'Moliere'
                f = electronScatteringFactorMoliere(Z(xx),Gmag(ii),E0*1e3);
        end
        F_hkl(ii) = F_hkl(ii) ...
            + f *exp(-2i*pi*(hkl(ii,1)*lattice(xx,1) + hkl(ii,2)*lattice(xx,2)...
            + hkl(ii,3)*lattice(xx,3)));       
    end
end

SFMag = abs(F_hkl); % Structure factor magnitude

%% Compute diffraction vs thickness

% Sample tilt
thetaX = 0; % rad, x component of tilt
thetaY = 0; % rad, y component of tilt
% Other inputs
Vcryst = prod(cellDim); % Unit cell volume in Angstroms^3
alpha = 0; % Inverse absorption depth (Angstroms^-1)
UCmax = 250;
dCrystal = (1:1:UCmax).*cellDim(3);
urms = 0.0894; % 1D rms displacement perpendicular to Bragg planes (Angstroms)
B = 8*pi^2*urms^2; % Debye-Waller B Factor

% Compute diffracted intensities vs thickness
[I,~] = computeDiffKin(...
    dCrystal,... % thicknesses
    thetaX,thetaY,... % sample tilt angles
    Ghkl,SFMag,... % g vectors and structure factors of reflections
    lamElec,... % electron wavelength
    Vcryst,... % unit cell volume
    alpha,... % absorption length
    B); % Debye-Waller Factor

Inorm = I./max(I,[],2); % Normalized intensities vs thickness

[savefile,savepath] = uiputfile('*.mat');
save([savepath savefile],'UCmax','thetaX','thetaY',...
    'I','dCrystal','Ghkl','SFMag','lamElec','Vcryst','alpha','B','cellDim')

%% Plot selected orders vs thickness
useNorm = false;

% List of q values of peaks to study
h1 = [0.4902 0]; 
h2 = [0 0.4902];
qxyPeaks = [1 0; ...
    1 1; 2 0; ...
    2 1; 2 2; ...
    3 0; 3 1; ...
    3 2; 2 3]...
    *[h1; h2];
peakNames = {'200',...
    '220','400',...
    '420','440',...
    '600','620',...
    '640'};

nPeaks = size(qxyPeaks,1);
indPeaks = zeros(nPeaks,1);
for iPeak = 1:nPeaks
    [~,indPeaks(iPeak)] = ...
        min((qxyPeaks(iPeak,1)-Ghkl(:,1)).^2 ...
        + (qxyPeaks(iPeak,2)-Ghkl(:,2)).^2);
end

nOrdersToShow = numel(indPeaks);
plotColors = jet(nOrdersToShow).*0.75;
orderNames = cell(nOrdersToShow,1);

figure;
for iShow = 1:nOrdersToShow
    orderShow = indPeaks(iShow); 
    if useNorm
        Iplot = Inorm(orderShow,:);
    else
        Iplot = I(orderShow,:);
    end
    plot(dCrystal./10,Iplot,'Color',plotColors(iShow,:)); hold all;  
    orderNames{iShow} = num2str(hkl(orderShow,:));
end
set(gca,'ytick',[],'fontsize',13); 
xlabel('thickness (nm)'); 
if useNorm
    ylabel('Normalized intensity (a.u.)');
else
    ylabel('Intensity (a.u.)')
end
legend(orderNames); legend boxoff

%% Diffraction + 2D plots at single thickness

% Sample tilt
thetaX = 0; % rad, x component of tilt
thetaY = 0; % rad, y component of tilt
% Other inputs
Vcryst = prod(cellDim); % Unit cell volume in Angstroms^3
alpha = 0; % Inverse absorption depth (Angstroms^-1)
nUC = 50;
t = cellDim(3)*nUC; % Crystal thickness (Angstroms)
B = 0; % Debye-Waller B Factor

% Compute diffracted intensities vs thickness
[I,excitationError] = computeDiffKin(...
    t,... % thicknesses
    thetaX,thetaY,... % sample tilt angles
    Ghkl,SFMag,... % g vectors and structure factors of reflections
    lamElec,... % electron wavelength
    Vcryst,... % unit cell volume
    alpha,... % absorption length
    B); % Debye-Waller Factor

% diffracted intensities
% NOTE: ASSUMES 001 ZONE AXIS, NEED TO PROJECT ONTO EWALD SPHERE IN GENERAL
figure;
scatter(Ghkl(:,1),Ghkl(:,2),50,I,'filled','MarkerEdgeColor','k')
xlabel('q_{x} (Angstroms^{-1})')
ylabel('q_{y} (Angstroms^{-1})')
title(['Intensity: t = ' num2str(t*0.1,3) ' nm, KE = ' num2str(E0) ' keV'])
colormap(violetFire)
colorbar
Isort = sort(I(:));
caxis([0 Isort(end-1)])

% excitation error
figure;
scatter(Ghkl(:,1),Ghkl(:,2),50,excitationError,'filled','MarkerEdgeColor','k')
xlabel('q_{x} (Angstroms^{-1})')
ylabel('q_{y} (Angstroms^{-1})')
title(['Excitation error: t = ' num2str(t*0.1,3) ' nm, KE = ' num2str(E0) ' keV'])
colormap('jet')
caxis([-1 1].*max(abs(excitationError)))
colorbar