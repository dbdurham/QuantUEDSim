%% Convergence tests -- Gold
clear all; clc;close all;

%% Set up simulation

tic
disp('Setting up...')
sDiff = setupSimBW();
toc

%% Choose diffracted beams to study

% hklTest = [2 0 0;...
%     2 2 0;...
%     4 0 0;...
%     4 2 0;...
%     4 4 0;...
%     6 0 0;...
%     6 2 0];
hklTest = [2 0 0;...
    -2 0 0;...
    2 2 0;...
    2 -2 0;...
    -2 2 0;...
    -2 -2 0;...
    4 0 0;...
    -4 0 0;...
    4 2 0;...
    -4 2 0;...
    4 -2 0;...
    -4 -2 0;...
    2 4 0;...
    -2 4 0;...
    2 -4 0;...
    -2 -4 0;...
    4 4 0;...
    4 -4 0;...
    -4 4 0;...
    -4 -4 0;...
    6 0 0;...
    -6 0 0;...
    6 2 0;...
    6 -2 0;...
    -6 2 0;...
    -6 -2 0;...
    -6 0 0;...
    2 6 0;...
    2 -6 0;...
    -2 6 0;...
    -2 -6 0;...
    ];
nPeaks = size(hklTest,1);

peakNames = cell(nPeaks,1);
for iPeak = 1:nPeaks
    peakNames{iPeak} = strrep(num2str(hklTest(iPeak,:)),' ','');
end

GhklTest = computeScatteringVectors(hklTest,sDiff.Gvec);

%% Perform simulations for varying params, extract ints vs thickness

% paramToTest = 'GxyThresh';
% paramRange = [2 2.5 3 3.5 4 4.5 5];
% sDiff.sThresh = 0.4/sDiff.cellDim(3);

% paramToTest = 'sThresh';
% paramRange = [0.1 0.15 0.2 0.3 0.4 0.5 0.6]./sDiff.cellDim(3);
% sDiff.GxyThresh = 4.5;

nTests = length(paramRange);

% Simulation parameters
nUCs = 40; % number of sim cells
theta_x = 0;
theta_y = 0.16; %rad

IArray = zeros(nPeaks,nUCs,nTests);
I0Array = zeros(nUCs,nTests);
tArray = 0.1*sDiff.cellDim(3)*(1:nUCs);

for iTest = 1:nTests
    % Update params
    switch paramToTest
        case 'GxyThresh'
            sDiff.GxyThresh = paramRange(iTest);
        case 'sThresh'
            sDiff.sThresh = paramRange(iTest);
    end
    
    rng(0) % fix random seed for comparison
    
    % Perform simulation
    tic
    disp(['Computing test # ' num2str(iTest)])
    [IArraySel,~,hklSel,GhklSel] = calcIntsBW(theta_x,theta_y,nUCs,...
        sDiff);
    % Extract zero beam intensity
    indZero = find(0 == hklSel(:,1) ...
            & 0 == hklSel(:,2) ...
            & 0 == hklSel(:,3));
    I0Array(:,iTest) = IArraySel(indZero,:);
    % Project diffracted ints onto DP
    NDP = [64 64];
    pixelSize = 2*sDiff.cellDim(1:2)./NDP;
    [qxa,qya] = makeFourierCoords(NDP,pixelSize);
    IDiff = projectIntsToDP(IArraySel,GhklSel,qxa,qya);
    % Extract diffracted ints being tested
    IArray(:,:,iTest) = extractIntsFromDP(IDiff,qxa,qya,GhklTest);
    toc

end

%% Plot intensity vs thickness for each peak as function of parameter

showIvtVsParam(IArray,tArray,peakNames,paramToTest,paramRange)

%% View the last generated diffraction pattern stack

StackViewerDiff(fftshift(fftshift(IDiff,1),2),tArray);

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

testToShow = 1;

showPctErrorPerPeak(IArray,tArray,testToShow,...
    peakNames,paramToTest,paramRange)

