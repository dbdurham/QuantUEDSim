%% Convergence tests -- Gold

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

% paramToTest = 'potBound';
% paramRange = [1 2 3 4 5];
% 
paramToTest = 'imageSizeCell';
paramRange = [64 128 256 512];

nTests = length(paramRange);

% Simulation parameters
nUCs = 100; % number of sim cells
theta_x = 0;
theta_y = 0; %rad

coefs = [nUCs,theta_x,theta_y];

IArray = zeros(nPeaks,nUCs,nTests);
I0Array = zeros(nUCs,nTests);

useGPU = false;

options = struct;
options.downSampFacCell = 1;
options.E0 = 4000e3;

for iTest = 1:nTests
    % Setup simulation
    
    options.(paramToTest) = paramRange(iTest);
    
    % Perform simulation
    sDiff = setupSimMS(options);
    if useGPU
        expPot = gpuArray(single(sDiff.expPot));
        tic
        [EWlast,EWstore,sDiff] = calcDiffMSGPU(sDiff,coefs,expPot);
        IDiff = double(abs(gather(EWstore)).^2);
        toc
    else
        [EWlast,EWstore,sDiff] = calcDiffMSCPU(sDiff,coefs);
        IDiff = double(abs(EWstore).^2);
    end
    I0Array(:,iTest) = IDiff(1,1,:);
    GhklTest = computeScatteringVectors(hklTest,sDiff.Gvec);
    IArray(:,:,iTest) = extractIntsFromDP(IDiff,...
        sDiff.qxaStore,sDiff.qyaStore,GhklTest);

end

%% Plot intensity vs thickness for each peak as function of parameter

tArray = 0.1*sDiff.cellDim(3)*(1:nUCs);
showIvtVsParam(IArray,tArray,peakNames,paramToTest,paramRange)

%% View the last generated diffraction pattern stack

StackViewerDiff(fftshift(fftshift(IDiff,1),2),tArray);

%% Compute and plot DP errors vs thickness for each parameter result

MPctE = computeMPctEStack(IArray);
R = computeRStack(IArray);
showErrorsVsThickness(tArray,R,MPctE,paramToTest,paramRange);

%% Compute and plot R for thickness bands

tBands = [0 5; 10 15; 15 20; 20 30; 30 40];
% tBands = [10 15];

RBands = computeRBands(IArray,tArray,tBands);

showRBands(RBands,tBands,paramRange,paramToTest);

%% Plot % error in each peak for a certain setting

testToShow = 2;

showPctErrorPerPeak(IArray,tArray,testToShow,...
    peakNames,paramToTest,paramRange)

