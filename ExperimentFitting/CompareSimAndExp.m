%% Show multislice diffraction patterns

% Function to generate two variables:
% hklExp: N x 3 array of hkl values for measured peaks
% GhklExp: N x 3 array of reciprocal space vectors
% intsExp: N x 1 array of measured intensities
[filename, pathname] = uigetfile('*.mat','Load IhklExp');
load([pathname filename]) 

%% Load stack of simulated patterns
[filename, pathname] = uigetfile('*.mat','Load DP library');
load([pathname filename])

Ilib1 = Ilib;

tArray = (1:nUC).*0.1*sDiff.cellDim(3);

% View pattern stack for max tilt spread
StackViewerDiff(fftshift(fftshift(Ilib1(:,:,:,end,iEnd),1),2),tArray)

%% Generate R maps between simulated and experimental patterns

Ny = size(Ilib1,1);
Nx = size(Ilib1,2);
nUC = size(Ilib1,3);
nTheta = size(Ilib1,4);

Rarray = zeros(nUC,nTheta);
IArrayExp = repmat(intsExp,[1,nUC]);
for iTheta = 1:nTheta
    IArraySim1 = extractIntsFromDP(Ilib1(:,:,:,iTheta,iEnd),...
        sDiff.qxaStore,sDiff.qyaStore,GhklExp);
    Rstack = computeRStack(cat(3,...
        IArraySim1,IArrayExp));
    Rarray(:,iTheta) = Rstack(:,1);
end

cmap = violetFire.^0.5;

figure;
imagesc(sigmaThetaSamp([1 end])*1e3,...
    tArray([1 end]),...
    Rarray);
xlabel('\sigma_{\theta} (mrad)')
ylabel('Thickness (nm)')
title('R_{exp - sim} (%)')
colormap(cmap)
colorbar()
caxis([0 60])
set(gca,'ydir','normal')

[~,indMin] = min(Rarray(:));
[subMin1,subMin2] = ind2sub(size(Rarray),indMin);

sigmaThetaMin = sigmaThetaSamp(subMin2)*1e3;
tMin = tArray(subMin1);

hold on
plot(sigmaThetaMin,tMin,...
    'o','Color',[0 0.5 0],...
    'MarkerSize',12,'LineWidth',1.5)

disp(['Best-fit thickness (nm): ' num2str(tMin,3) ])
disp(['Best-fit RMS tilt spread (mrad): ' num2str(sigmaThetaMin,3)])
disp(['Best residual (%): ' num2str(Rarray(subMin1,subMin2))])

%% Show simulated intensities
iUC = subMin1;
iTheta = subMin2;

% Extract sim intensities
IArraySim1 = extractIntsFromDP(Ilib1(:,:,:,iTheta,iEnd),...
        sDiff.qxaStore,sDiff.qyaStore,GhklExp);
intsSim = IArraySim1(:,iUC);

% Plot intensity vs thickness for the max iteration
I0ArraySim1 = extractIntsFromDP(Ilib1(:,:,:,iTheta,iEnd),...
        sDiff.qxaStore,sDiff.qyaStore,[0 0 0]);
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

%% Load a second simulated library to compare against the first

[filename, pathname] = uigetfile('*.mat','Load DP library');
load([pathname filename])
Ilib2 = Ilib;
tArray2 = (1:nUC).*0.1*sDiff.cellDim(3);

% View pattern stack for max tilt spread
StackViewerDiff(fftshift(fftshift(Ilib2(:,:,:,end,iEnd),1),2),tArray2)

%% Plot intensity vs thickness for max iteration

IArraySim2 = extractIntsFromDP(Ilib2(:,:,:,iTheta,iEnd),...
        sDiff.qxaStore,sDiff.qyaStore,GhklExp);
I0ArraySim2 = extractIntsFromDP(Ilib2(:,:,:,iTheta,iEnd),...
        sDiff.qxaStore,sDiff.qyaStore,[0 0 0]);
showIvt(IArraySim2,I0ArraySim2,tArray,peakNames);

%% Generate R map comparing the two simulated datasets

Rarray = zeros(nUC,nTheta);
for iTheta = 1:nTheta
    IArraySim2 = extractIntsFromDP(Ilib2(:,:,:,iTheta,iEnd),...
        sDiff.qxaStore,sDiff.qyaStore,GhklExp);
    IArraySim1 = extractIntsFromDP(Ilib1(:,:,:,iTheta,iEnd),...
        sDiff.qxaStore,sDiff.qyaStore,GhklExp);
    Rstack = computeRStack(cat(3,...
        IArraySim1,IArraySim2));
    Rarray(:,iTheta) = Rstack(:,1);
end

cmap = violetFire.^0.5;

figure;
imagesc(sigmaThetaSamp([1 end])*1e3,...
    tArray([1 end]),...
    Rarray);
xlabel('\sigma_{\theta} (mrad)')
ylabel('Thickness (nm)')
title('R_{sim2 - sim1} (%)')
colormap(cmap)
colorbar()
caxis([0 60])
set(gca,'ydir','normal')

%% UPDATE OR REMOVE THE CODE BELOW 


















%% Load kinematical data

load('KinDiff_Gold_750keV_160mradTilt_10Iter.mat')

u = 0.0894;

h1 = [0.4902 0]; 
h2 = [0 0.4902];
qxyPeaks = [1 0; ...
    1 1; 2 0; ...
    2 1; 2 2; ...
    3 0; 3 1]...
    *[h1; h2];
peakNames = {'200',...
    '220','400',...
    '420','440',...
    '600','620'};

nPeaks = size(qxyPeaks,1);
indPeaks = zeros(nPeaks,1);
for iPeak = 1:nPeaks
    [~,indPeaks(iPeak)] = ...
        min((qxyPeaks(iPeak,1)-Ghkl(:,1)).^2 ...
        + (qxyPeaks(iPeak,2)-Ghkl(:,2)).^2);
end

IarrayKin = Ilib(indPeaks,:,:,iEnd);
% Apply Debye-Waller Factor
q2 = qxyPeaks(:,1).^2 + qxyPeaks(:,2).^2;
IarrayKin = IarrayKin.*exp(-4*pi^2*u^2*q2);

%% R factor map kinematical

RarrayKin = zeros(nUC,nTheta);
for iUC = 1:nUC
    for iTheta = 1:nTheta
        intsSim = IarrayKin(:,iUC,iTheta);
        % Adjust scale for best fit
        scaleFac = intsExp.^gam\intsSim.^gam;
        intsSim = (intsSim.^gam./scaleFac).^(1/gam);
        RarrayKin(iUC,iTheta) = ...
            sum(abs(intsSim.^gam - intsExp.^gam))./sum(intsExp.^gam);
    end
end

% Generate custom colormap for error mapping
cpts = [1         1         1        0;
        218/255   165/255   32/255   0.2;
        0         0         0        1.0];
cmap = generateGradColormap(cpts,1024);

figure;
imagesc(sigmaThetaSamp([1 end])*1e3,...
    [1 nUC].*0.1*sDiff.cellDim(3),...
    100*RarrayKin);
xlabel('\sigma_{\theta} (mrad)')
ylabel('Thickness (nm)')
title('R_{exp - kin} (%)')
colormap(cmap)
colorbar()
caxis([0 60])

% hold on
% plot(77.5,12.24,'o','Color',[0 0.5 0],...
%     'MarkerSize',12,'LineWidth',1.5)

%% Scatter plot kinematical + multislice

iUC = subMin1;
iTheta = subMin2;

% Multislice
intsSim = zeros(size(intsExp));
% Extract sim intensities
for iPeak = 1:size(hklExp,1)
    indSim = find(hklExp(iPeak,1) == hklSim(:,1) & ...
        hklExp(iPeak,2) == hklSim(:,2) & ...
        hklExp(iPeak,3) == hklSim(:,3));
    [subSim1,subSim2] = ind2sub([Ny,Nx],indSim);
    intsSim(iPeak) = Ilib(subSim1,subSim2,iUC,iTheta,8);
end
% Adjust scale for best fit
scaleFac = intsExp.^gam\intsSim.^gam;
intsSim = (intsSim.^gam./scaleFac).^(1/gam);

% Kinematical
intsKin = IarrayKin(:,iUC,iTheta);
scaleFac = intsExp.^gam\intsKin.^gam;
intsKin = (intsKin.^gam./scaleFac).^(1/gam);

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
plot(intsExp.^gam,intsSim.^gam,'b.','MarkerSize',12)
hold on
plot(intsExp.^gam,intsKin.^gam,'r.','MarkerSize',12)
plot([0 3e3].^gam,[0 3e3].^gam,'k--')
xlabel(['I_{exp}^{' num2str(gam) '}'])
ylabel(['I_{sim}^{' num2str(gam) '}'])
xlim([0 3e3].^gam);
ylim([0 3e3].^gam);
axis square
legend('MS','Kin')