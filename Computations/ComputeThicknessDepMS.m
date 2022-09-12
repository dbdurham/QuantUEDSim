%% Compute thickness dependent diffraction peaks + wave functions

%% 1. Set up simulation cell + parameters

% Can set parameters here through the options struct, 
% or in the setup function itself
options = struct;
options.E0 = 750e3; % Electron energy (eV)
options.uRMS = 0.0894; % 1D rms atomic displacement
options.qRange = 4; % q range to store in *output* (A^-1)
options.cellMult = 1;
options.downSampFacCell = 1;
sDiff = setupSimMS(options); 

%% 2a. Compute diffraction patterns
sDiff.storeRealSpaceCellEWs = false;
nUCs = 100; % Number of simulation cells to propagate through
theta1 = 0.0; % Incident angle of the e-beam (rad)
theta2 = 0.0;
IDiff = calcDiffMS(theta1,theta2,nUCs,sDiff);

tArray = 0.1*sDiff.cellDim(3)*(1:nUCs);
StackViewerDiff(fftshift(fftshift(IDiff,1),2),tArray)

%% 2b. Visualize diffraction peaks

hklTest = [2 0 0;...
    2 2 0;...
    4 0 0;...
    4 2 0;...
    4 4 0;...
    6 0 0;...
    6 2 0];

nPeaks = size(hklTest,1);

peakNames = cell(nPeaks,1);
for iPeak = 1:nPeaks
    peakNames{iPeak} = strrep(num2str(hklTest(iPeak,:)),' ','');
end

GhklTest = computeScatteringVectors(hklTest,sDiff.Gvec);

IArray = extractIntsFromDP(IDiff,...
    sDiff.qxaStore,sDiff.qyaStore,GhklTest);
I0Array = extractIntsFromDP(IDiff,...
    sDiff.qxaStore,sDiff.qyaStore,[0 0 0]);

showIvt(IArray,I0Array,tArray,peakNames);

%% 2c. Compute and overlay kinematical diff signals
optionsKin.E0 = options.E0;
optionsKin.uRMS = options.uRMS;
sDiffKin = setupSimKin(optionsKin);
IDiffKin = calcDiffKin(theta1,theta2,nUCs,sDiffKin);

IArrayKin = extractIntsFromDP(IDiffKin,...
    sDiffKin.qxaStore,sDiffKin.qyaStore,GhklTest);

showIvtPlusKin(IArray,I0Array,...
    IArrayKin,tArray,peakNames);

%% 2d. Save computations

save(['DiffAnalysis_' num2str(sDiff.E0/1e3) 'kV.mat'],...
    'tArray','Rarray','IArray','I0Array','IArrayKin',...
    'sDiff','sDiffKin')

%% 3a. Compute real space electron waves
sDiff.storeRealSpaceCellEWs = true;
nUCs = 100; % Number of simulation cells to propagate through
theta1 = 0.0; % Incident angle of the e-beam (rad)
theta2 = 0.0;
coefs = [nUCs,theta1,theta2];
[~,EWImage,sDiff] = calcDiffMSCPU(sDiff,coefs);

tArray = 0.1*sDiff.cellDim(3)*(1:nUCs); % Thickness (nm)
StackViewerEWs(abs(EWImage),tArray)

%% 3b. Visualize real space waves

cpts = [0 0 0 0;
    0 0.5 0 0.24;
    220/255 202/255 152/255 0.8;
    0.95 0.95 0.8 1];
N = 1024;
cmapMag = generateGradColormap(cpts,N);

figure('Position',[100 100 350 400]);
imData = squeeze(abs(EWImage(64,:,:)))';
imagesc([0 0.1*sDiff.cellDim(1)],[tArray(1) tArray(end)],...
    imData)
xlabel('x (nm)')
ylabel('z (nm)')
title('|\psi|')
colorbar
colormap(cmapMag)
caxis([0 6.2])

figure('Position',[100 100 350 400]);
imagesc([0 0.1*sDiff.cellDim(1)],[tArray(1) tArray(end)],...
    squeeze(angle(EWImage(64,:,:)))')
xlabel('x (nm)')
ylabel('z (nm)')
title('\phi')
colorbar
colormap(hsv(1024))
