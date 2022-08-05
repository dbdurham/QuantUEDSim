%% Compute thickness dependent diffraction

simType = 'Multislice'; %Simulation method: 
% Kinematical, Bloch Waves, or Multislice
theta1 = 0; % rad, x component of tilt
theta2 = 0; % rad, y component of tilt
nUC = 100; % Number of unit cells to simulate

% Setup simulation
sDiff = setupSim(simType);
% Compute diffraction patterns
IDiff = calcDiff(theta1,theta2,nUC,sDiff);

tArray = (1:nUC)*sDiff.cellDim(3); % Angstroms
StackViewerDiff(fftshift(fftshift(IDiff,1),2),tArray.*0.1)

%% Extract and plot diffraction peak intensities vs thickness

% Diffraction orders to show
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

I0Array = extractIntsFromDP(IDiff,...
    sDiff.qxaStore,sDiff.qyaStore,[0 0 0]);
IArray = extractIntsFromDP(IDiff,...
    sDiff.qxaStore,sDiff.qyaStore,GhklTest);

% Plot diffracted intensities vs thickness
showIvt(IArray,I0Array,tArray,peakNames);