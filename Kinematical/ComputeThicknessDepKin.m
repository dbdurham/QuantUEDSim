%% Thickness dependent kinematical calculation

% Set up Bloch Wave simulations
sDiff = setupSimKin;

% Inputs for calculation
theta1 = 0; % rad, x component of tilt
theta2 = 0; % rad, y component of tilt
nUC = 100; % Number of unit cells to simulate

% Calculate diffraction pattern stack
IDiff = calcDiffKin(theta1,theta2,nUC,sDiff);

tArray = (1:nUC)*sDiff.cellDim(3)*0.1; % nm
% View the diffraction patterns
StackViewerDiff(fftshift(fftshift(IDiff,1),2),tArray)

%% Extract intensities of interest

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

IArray = extractIntsFromDP(IDiff,sDiff.qxaStore,sDiff.qyaStore,GhklTest);
I0Array = ones(1,size(IArray,2));

% Plot diffracted intensities vs thickness
showIvt(IArray,I0Array,tArray,peakNames)