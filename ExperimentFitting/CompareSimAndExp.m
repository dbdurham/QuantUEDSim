%% Compare dynamical scattering library to experimental static pattern
% Figure 5c in D.B. Durham et al, ArXiv Preprint 2022


%% Load experimental intensities
% Function to generate two variables:
% hklExp: N x 3 array of hkl values for measured peaks
% GhklExp: N x 3 array of reciprocal space vectors
% intsExp: N x 1 array of measured intensities
[filename, pathname] = uigetfile('*.mat','Load IhklExp');
load([pathname filename]) 

%% Load stack of simulated patterns
[filename, pathname] = uigetfile('*.mat','Load DP library');
load([pathname filename])

sDiff1 = sDiff;
Ilib1 = Ilib;

tArray = (1:nUC).*0.1*sDiff1.cellDim(3);

% View pattern stack for max tilt spread
StackViewerDiff(fftshift(fftshift(Ilib1(:,:,:,end,end),1),2),tArray)

%% Generate R maps between simulated and experimental patterns

Ny = size(Ilib1,1);
Nx = size(Ilib1,2);
nUC = size(Ilib1,3);
nTheta = size(Ilib1,4);

% applyDWF = false;
% uRMS = 0.0894;
% [~,DWFInt] = computeDWF(uRMS,1,GmagExp);

Rarray = zeros(nUC,nTheta);
IArrayExp = repmat(intsExp,[1,nUC]);
for iTheta = 1:nTheta
    IArraySim1 = extractIntsFromDP(Ilib1(:,:,:,iTheta,end),...
        sDiff1.qxaStore,sDiff1.qyaStore,GhklExp);
    if applyDWF
        IArraySim1 = IArraySim1.*DWFInt;
    end
    Rstack = computeRStack(cat(3,...
        IArraySim1,IArrayExp));
    Rarray(:,iTheta) = Rstack(:,1);
end

cmap = violetFire.^0.5;

figure('Position',[200 200 280 200]);
imagesc(sigmaThetaSamp([1 end])*1e3,...
    tArray([1 end]),...
    Rarray);
xlabel('\sigma_{\theta} (mrad)')
ylabel('Thickness (nm)')
title('R_{Exp - Sim} (%)')
colormap(cmap)
colorbar()
caxis([0 50])
set(gca,'ydir','normal')

[~,indMin] = min(Rarray(:));
[subMin1,subMin2] = ind2sub(size(Rarray),indMin);

sigmaThetaMin = sigmaThetaSamp(subMin2)*1e3;
tMin = tArray(subMin1);

hold on
plot(sigmaThetaMin,tMin,...
    'o','Color',[0 0.5 0],...
    'MarkerSize',8,'LineWidth',1)

disp(['Best-fit thickness (nm): ' num2str(tMin,3) ])
disp(['Best-fit RMS tilt spread (mrad): ' num2str(sigmaThetaMin,3)])
disp(['Best residual (%): ' num2str(Rarray(subMin1,subMin2))])

%% Show simulated intensities
iUC = subMin1;
iTheta = subMin2;

% Extract sim intensities
IArraySim1 = extractIntsFromDP(Ilib1(:,:,:,iTheta,end),...
        sDiff1.qxaStore,sDiff1.qyaStore,GhklExp);
intsSim = IArraySim1(:,iUC);

% Plot intensity vs thickness for the max iteration
I0ArraySim1 = extractIntsFromDP(Ilib1(:,:,:,iTheta,end),...
        sDiff1.qxaStore,sDiff1.qyaStore,[0 0 0]);
showIvt(IArraySim1,I0ArraySim1,tArray,peakNames);

% Adjust scale for best fit
gam = 0.5;
scaleFac = intsExp.^gam\intsSim.^gam;
intsSim = (intsSim.^gam./scaleFac).^(1/gam);

% Show the avg off intensities
figure;
scatter(hklExp(:,1),hklExp(:,2),intsSim.^0.5,...
    intsSim.^0.5,'Filled','MarkerEdgeColor','k');
colormap(parula(1024))
colorbar()
set(gca,'ydir','reverse')
axis equal off
title('Simulated static pattern I^{-1/2}')

% Scatter plot of simulated vs measured intensities
figure('Position',[200 200 500 400]);
scatter(intsExp.^gam,intsSim.^gam,200,'b.')
hold on
plot([0 2e3].^gam,[0 2e3].^gam,'k--')
xlabel(['I_{exp}^{' num2str(gam) '}'])
ylabel(['I_{sim}^{' num2str(gam) '}'])
xlim([0 2e3].^gam);
ylim([0 2e3].^gam);
axis square

%% Load kinematical data

[filename, pathname] = uigetfile('*.mat','Load Kinematical data');
load([pathname filename])

sDiffKin = sDiff;
IlibKin = Ilib;

% Apply Debye-Waller Factor
GmagStore = sqrt(sDiffKin.qxaStore.^2 + sDiffKin.qyaStore.^2);
[~,DWFInt,~] = computeDWF(sDiff1.uRMS,1,GmagStore);
IlibKin = IlibKin.*DWFInt;

nPeaks = size(GhklExp,1);
IArrayKin = zeros(nPeaks,nUC,nTheta);
for iTheta = 1:nTheta
    IArrayKin(:,:,iTheta) = extractIntsFromDP(IlibKin(:,:,:,iTheta,end),...
        sDiffKin.qxaStore,sDiffKin.qyaStore,GhklExp);
end

%% R factor map kinematical

RarrayKin = zeros(nUC,nTheta);
for iUC = 1:nUC
    for iTheta = 1:nTheta
        intsSim = IArrayKin(:,iUC,iTheta);
        % Adjust scale for best fit
        scaleFac = intsExp.^gam\intsSim.^gam;
        intsSim = (intsSim.^gam./scaleFac).^(1/gam);
        RarrayKin(iUC,iTheta) = ...
            sum(abs(intsSim.^gam - intsExp.^gam))./sum(intsExp.^gam);
    end
end

cmap = violetFire.^0.5;

figure('Position',[200 200 280 200]);
imagesc(sigmaThetaSamp([1 end])*1e3,...
    [1 nUC].*0.1*sDiff.cellDim(3),...
    100*RarrayKin);
xlabel('\sigma_{\theta} (mrad)')
ylabel('Thickness (nm)')
title('R_{exp - kin} (%)')
set(gca,'ydir','normal')
colormap(cmap)
colorbar()
caxis([0 50])

% hold on
% plot(77.5,12.24,'o','Color',[0 0.5 0],...
%     'MarkerSize',12,'LineWidth',1.5)