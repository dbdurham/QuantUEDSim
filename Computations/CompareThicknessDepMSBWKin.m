%% Compare thickness dependence: BW vs MS vs Kin
% as used in D.B. Durham et al, ArXiv Preprint 2022 

% Set up simulations
sDiffBW = setupSimBW;
sDiffMS = setupSimMS;

% Inputs for calculation
theta1 = 0; % rad, x component of tilt
theta2 = 0; % rad, y component of tilt
nUC = 100; % Number of unit cells to simulate

IDiffBW = calcDiffBW(-theta1,-theta2,nUC,sDiffBW);
IDiffMS = calcDiffMS(theta1,theta2,nUC,sDiffMS);

%% Display and compare

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

GhklTest = computeScatteringVectors(hklTest,sDiffBW.Gvec);

IArrayBW = extractIntsFromDP(IDiffBW,...
    sDiffBW.qxaStore,sDiffBW.qyaStore,GhklTest);
I0ArrayBW = extractIntsFromDP(IDiffBW,...
    sDiffBW.qxaStore,sDiffBW.qyaStore,[0 0 0]);
IArrayMS = extractIntsFromDP(IDiffMS,...
    sDiffMS.qxaStore,sDiffMS.qyaStore,GhklTest);
I0ArrayMS = extractIntsFromDP(IDiffMS,...
    sDiffMS.qxaStore,sDiffMS.qyaStore,[0 0 0]);

tArray = (1:nUC)*sDiffBW.cellDim(3)*0.1; %nm

showIvt(IArrayBW,I0ArrayBW,tArray,peakNames);
title('Bloch Waves')
showIvt(IArrayMS,I0ArrayMS,tArray,peakNames);
title('Multislice')

R = computeRStack(cat(3,IArrayBW,IArrayMS));
figure;
semilogy(tArray,R(:,1),'k-','LineWidth',1.5)
xlabel('Thickness (nm)')
ylabel('R_{BW-MS} (%)')

%% Plot BW intensities superimposed w/ Kinematical intensities

sDiffKin = setupSimKin;
IDiffKin = calcDiffKin(-theta1,-theta2,nUC,sDiffKin);
% Apply Debye-Waller Factor
GmagStore = sqrt(sDiffKin.qxaStore.^2 + sDiffKin.qyaStore.^2);
[~,DWFInt,~] = computeDWF(sDiffBW.uRMS,1,GmagStore);
IDiffKin = IDiffKin.*DWFInt;

IArrayKin = extractIntsFromDP(IDiffKin,...
    sDiffKin.qxaStore,sDiffKin.qyaStore,GhklTest);
I0ArrayKin = extractIntsFromDP(IDiffKin,...
    sDiffKin.qxaStore,sDiffKin.qyaStore,[0 0 0]);
showIvtPlusKin(IArrayBW,I0ArrayBW,...
    IArrayKin,tArray,peakNames);