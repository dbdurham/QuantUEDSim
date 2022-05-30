%% Convergence tests -- Gold
clear all; clc;close all;
%% Choose diffracted beams to study

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

%% Perform simulations for varying params, extract ints vs thickness

% paramToTest = 'cellMult';
% paramRange = [1 2 4];
% isWyckoffOpt = false;
% 
% paramToTest = 'potBound';
% paramRange = [1 2 3 4 5];
% isWyckoffOpt = false;
% 
paramToTest = 'imageSizeCell';
paramRange = [64 128 256 512];
isWyckoffOpt = false;

nTests = length(paramRange);

% Simulation parameters
nUCs = 40; % number of sim cells
theta_x = 0;
theta_y = 0; %rad

IArray = zeros(nPeaks,nUCs,nTests);
I0Array = zeros(nUCs,nTests);

useGPU = false;

for iTest = 1:nTests
    % Setup simulation
    options = struct;
    options.(paramToTest) = paramRange(iTest);
    
    rng(0) % fix random seed for comparison
    
    % Perform simulation
    tic
    sDiff = setupMultisliceSim(options);
    if useGPU
        [EWlast,EWstore,sDiff] = calcDiffGPU(sDiff,coefs);
        IDiff = double(abs(gather(EWstore)).^2);
    else
        [EWlast,EWstore,sDiff] = calcDiff(sDiff,coefs);
        IDiff = double(abs(EWstore).^2);
    end
    I0Array(:,iTest) = IDiff(1,1,:);
    GhklTest = computeScatteringVectors(hklTest,sDiff.Gvec);
    IArray(:,:,iTest) = extractIntsFromDP(IDiff,...
        sDiff.qxa,sDiff.qya,GhklTest);
    toc

end

%% Plot intensity vs thickness for each peak as function of parameter

tArray = 0.1*sDiff.cellDim(3)*(1:nUCs);
showIvtVsParam(IArray,tArray,peakNames,paramToTest,paramRange)

%% Compute and plot DP errors vs thickness for each parameter result

MPctE = computeMPctEStack(IArray);
R = computeRStack(IArray);
showErrorsVsThickness(tArray,R,MPctE,paramToTest,paramRange);

%% Compute and plot R for thickness bands

tBands = [0 5; 10 15; 15 20];
% tBands = [10 15];

RBands = computeRBands(IArray,tArray,tBands);

showRBands(RBands,tBands,paramRange,paramToTest);

%% Plot % error in each peak for a certain setting

testToShow = 2;

showPctErrorPerPeak(IArray,tArray,testToShow,...
    peakNames,paramToTest,paramRange)

