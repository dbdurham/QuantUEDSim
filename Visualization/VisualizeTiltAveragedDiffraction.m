%% Visualize tilt-averaged diffraction

%% Load tilt-averaged diffraction library

[loadfile,loadpath] = uigetfile('*.mat');
load([loadpath loadfile]);

%% Visualize resulting intensities and R factors

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

IArray = zeros(nPeaks,nUC,nIter);
I0Array = zeros(1,nUC,nIter);
iTheta = 16;
for iIter = 1:nIter
    IArray(:,:,iIter) = extractIntsFromDP(Ilib(:,:,:,iTheta,iIter),...
        sDiff.qxaStore,sDiff.qyaStore,GhklTest);
    I0Array(:,:,iIter) = extractIntsFromDP(Ilib(:,:,:,iTheta,iIter),...
        sDiff.qxaStore,sDiff.qyaStore,[0 0 0]);
end

%% Plot intensity vs thickness for each peak vs iteration

showIvtVsParam(IArray,tArray,peakNames,'Iterations',1:nIter)

%% Plot intensity vs thickness for the max iteration
showIvt(IArray(:,:,end),I0Array(:,:,end),tArray,peakNames);

%% Compute and plot R for thickness bands

tBands = [0 5; 10 15; 15 20; 25 30; 35 40];
% tBands = [10 15];

RBands = computeRBands(IArray,tArray,tBands);

showRBands(RBands,tBands,1:nIter,'Iterations');
