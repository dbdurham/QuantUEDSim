%% Bloch Waves diffraction calculation example

% Set up Bloch Wave simulations
sDiff = setupSimBW;

% Inputs for calculation
theta1 = 0; % rad, x component of tilt
theta2 = 0; % rad, y component of tilt
nUC = 100; % Number of unit cells to simulate
UThresh = 5e-3; % Lower threshold on scattering potential of beams to include
sThresh = 5e-2; % Upper threshold on excitation error of beams to include

% Calculate diffracted intensities for selected thicknesses
[psi_G_array,hklSel,isSel] = calcIntsBW(theta1,theta2,nUC,...
    UThresh,sThresh,sDiff);
tArray = (1:nUC)*sDiff.cellDim(3); % Angstroms

% Plot diffracted intensities vs thickness
figure;
plot(tArray.*0.1,abs(psi_G_array).^2,'-')
xlabel('Distance (nm)')
ylabel('Beam amplitude')

% Map diffracted intensities onto 2D diffraction pattern 
[IDiff,qxa,qya] = calcDiffBW(psi_G_array,isSel,sDiff);
% View the diffraction patterns
StackViewerDiff(fftshift(fftshift(IDiff,1),2),tArray.*0.1)