%% Compare tilt-averaged BW, MS

%% Load stack of BW simulated patterns
[filename, pathname] = uigetfile('*.mat','Load DP library');
load([pathname filename])

sDiff2 = sDiff;
Ilib2 = Ilib;
iEnd2 = size(Ilib2,5);

%% Load stack of MS simulated patterns

[filename, pathname] = uigetfile('*.mat','Load DP library');
load([pathname filename])

sDiff1 = sDiff;
Ilib1 = Ilib;
iEnd = size(Ilib1,5);

tArray = (1:nUC).*0.1*sDiff1.cellDim(3);

%% Identify peaks to compare

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

%% Generate R map comparing MS, BW

Rarray = zeros(nUC,nTheta);
for iTheta = 1:nTheta
    IArraySim2 = extractIntsFromDP(Ilib2(:,:,:,iTheta,iEnd2),...
        sDiff2.qxaStore,sDiff2.qyaStore,GhklTest);
    IArraySim1 = extractIntsFromDP(Ilib1(:,:,:,iTheta,iEnd),...
        sDiff1.qxaStore,sDiff1.qyaStore,GhklTest);
    Rstack = computeRStack(cat(3,...
        IArraySim1,IArraySim2));
    Rarray(:,iTheta) = Rstack(:,1);
end

cmap = violetFire.^0.5;

figure;
imagesc(sigmaThetaSamp([1 end])*1e3,...
    tArray([1 end]),...
    Rarray);
hold on
[X,Y] = meshgrid(sigmaThetaSamp.*1e3,tArray);
%contour(X,Y,Rarray,[3,3],'w--','LineWidth',1)
%contour(X,Y,Rarray,[10,10],'w--','LineWidth',1)
xlabel('\sigma_{\theta} (mrad)')
ylabel('Thickness (nm)')
title('R_{BW - MS} (%)')
colormap(cmap)
colorbar()
caxis([0 40])
set(gca,'ydir','normal')

%% Plot tilt-averaged peaks

iTheta = 40;
IArraySim2 = extractIntsFromDP(Ilib2(:,:,:,iTheta,iEnd2),...
        sDiff2.qxaStore,sDiff2.qyaStore,GhklTest);
I0ArraySim2 = extractIntsFromDP(Ilib2(:,:,:,iTheta,iEnd2),...
        sDiff2.qxaStore,sDiff2.qyaStore,[0 0 0]);
IArraySim1 = extractIntsFromDP(Ilib1(:,:,:,iTheta,iEnd),...
    sDiff1.qxaStore,sDiff1.qyaStore,GhklTest);

showIvtPlusKin(IArraySim2,I0ArraySim2,IArraySim1,tArray,peakNames);