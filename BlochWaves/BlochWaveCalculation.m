%% Bloch Waves diffraction calculation example

% Set up Bloch Wave simulations
sDiff = setupSimBW;

% Inputs for calculation
theta1 = 0; % rad, x component of tilt
theta2 = 0; % rad, y component of tilt
nUC = 100; % Number of unit cells to simulate

% Calculate diffracted intensities for selected thicknesses
[IArray,psi_G_array,hklSel,GhklSel] = calcIntsBW(theta1,theta2,nUC,...
    sDiff);
tArray = (1:nUC)*sDiff.cellDim(3); % Angstroms

% Plot diffracted intensities vs thickness
figure;
plot(tArray.*0.1,IArray,'-')
xlabel('Distance (nm)')
ylabel('Beam intensity')

% Map diffracted intensities onto 2D diffraction pattern
NDP = [32 32];
pixelSize = sDiff.cellDim(1:2)./NDP;
[qxa,qya] = makeFourierCoords(NDP,pixelSize);
IDiff = projectIntsToDP(IArray,GhklSel,qxa,qya);
% View the diffraction patterns
StackViewerDiff(fftshift(fftshift(IDiff,1),2),tArray.*0.1)

% % Test intensity extraction
% IArrayRet = extractIntsFromDP(IDiff,qxa,qya,GhklSel);