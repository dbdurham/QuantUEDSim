%% Convergence tests -- Gold
clear all; clc;close all;
%% Choose diffracted beams to study

% list q values of peaks to study
% % define some basis q vectors
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

%% Perform simulations for multiple parameters, 
% extract diffracted intensities vs thickness
% paramToTest = 'cellMult';
% paramRange = [1 2 4];
% isWyckoffOpt = false;
% 
% paramToTest = 'potBound';
% paramRange = [1 2 3 4 5];
% isWyckoffOpt = false;
% 
paramToTest = 'imageSizeCell';
paramRange = [64 128 256 512 1024];
isWyckoffOpt = false;

nTests = length(paramRange);

% Simulation parameters
nUCs = 150; % number of sim cells
theta_x = 0;
theta_y = 0; %rad
coefs = [nUCs,0,theta_x,theta_y];

intValsArray = zeros(nPeaks,nUCs,nTests);
I0Array = zeros(nUCs,nTests);

useGPU = false;

for iTest = 1:nTests
    % Setup simulation
    if isWyckoffOpt
        options = struct;
        wyckoffOptions = struct;
        wyckoffOptions.(paramToTest) = paramRange(iTest);
        sRefine = setupMultisliceSim(options,wyckoffOptions);
    else
        options = struct;
        options.(paramToTest) = paramRange(iTest);
        sRefine = setupMultisliceSim(options);
    end
    
    rng(0) % fix random seed for comparison
    
    % Perform multislice
    tic
    if useGPU
        [EWlast,EWstore,sRefine] = calcDiffGPU(sRefine,coefs);
        Istore = double(abs(gather(EWstore)).^2);
    else
        [EWlast,EWstore,sRefine] = calcDiff(sRefine,coefs);
        Istore = double(abs(EWstore).^2);
    end
    
    toc
    
    % Extract peak intensities vs thickness   
    indPeaks = zeros(nPeaks,1);
    for iPeak = 1:nPeaks
        [~,indPeaks(iPeak)] = ...
            min((qxyPeaks(iPeak,2)-sRefine.qxaStore(:)).^2 ...
            + (qxyPeaks(iPeak,1)-sRefine.qyaStore(:)).^2);
    end
    
    I0Array = zeros(1,nUCs);
    for iUC = 1:nUCs
        I = Istore(:,:,iUC);
        I0 = I(1,1);
        I0Array(iUC,iTest) = I0;
        % normalize peaks by total *diffracted* intensity
        % I = I./(sum(I(:))-I0); 
        intValsArray(:,iUC,iTest) = I(indPeaks);
    end
end

%% Plot intensity vs thickness for each peak as function of parameter
if strcmp(paramToTest,'imageSizeCell')
    legendStr = arrayfun(@(x) ['N = ' num2str(2*x) ' px'],paramRange,...
        'UniformOutput',false);
else
    legendStr = arrayfun(@(x) [paramToTest ' = ' num2str(x)],paramRange,...
        'UniformOutput',false);
end
figure
indsToPlot = 1:7;
nPlots = numel(indsToPlot);
for iPlot = 1:nPlots
    subplot(2,ceil(nPlots/2),iPlot)
    indPeak = indsToPlot(iPlot);
    plot(0.1*sRefine.cellDim(3)*(1:nUCs),...
        squeeze(intValsArray(indPeak,:,:)),...
        '-o')
    xlabel('Thickness (nm)')
    ylabel('Fraction of diffracted intensity')
    title(peakNames{indPeak})
    if iPlot == nPlots
       lg= legend(legendStr);
    end
end
% 
% set(gcf,'color','white','position',[50 50 850 400]);
% set(lg,'position',[0.85 0.35 0.1 0.1]);

%% Plot DP errors vs thickness for each parameter result

if strcmp(paramToTest,'imageSizeCell')
    legendStr = arrayfun(@(x) ['N = ' num2str(2*x) ' px'],paramRange,...
        'UniformOutput',false);
else
    legendStr = arrayfun(@(x) [paramToTest ' = ' num2str(x)],paramRange,...
        'UniformOutput',false);
end
tArray = 0.1*sRefine.cellDim(3)*(1:nUCs);

MSE = squeeze(mean(abs(diff(intValsArray,1,3)).^2,1));
MPctE = 100*squeeze(mean(abs(diff(intValsArray,1,3))...
    ./intValsArray(:,:,1:end-1),1));
R = 100*squeeze(sum(abs(diff(intValsArray,1,3)),1) ...
    ./ sum(intValsArray(:,:,1:end-1),1));

figure;
semilogy(tArray,R,'LineWidth',1.5)
xlabel('Thickness (nm)')
ylabel('R (%)')
legend(legendStr{1:end-1})

figure;
semilogy(tArray,MPctE,'LineWidth',1.5)
xlabel('Thickness (nm)')
ylabel('Mean % Error')
legend(legendStr{1:end-1})

%% Plot R for thickness bands

% tRanges = [0 5; 10 15; 15 20; 35 40; 55 60; 75 80; 95 100];
tRanges = [10 15];
tMean = mean(tRanges,2);
nRanges = size(tRanges,1);
colorList = 0.8*jet(nRanges);

Rend = 100*squeeze(...
    sum(abs(intValsArray(:,:,1:end) - intValsArray(:,:,end)),1) ...
    ./ sum(intValsArray(:,:,end),1));

figure;
for iRange = 1:nRanges
    indLow = find(tArray > tRanges(iRange,1),1,'first');
    indHigh = find(tArray < tRanges(iRange,2),1,'last');
    if strcmp(paramToTest,'imageSizeCell')
        xData = paramRange*sRefine.cellMult;
    else
        xData = paramRange;
    end   
    if nRanges == 1
        plot(xData,mean(Rend(indLow:indHigh,:),1),'k.-',...
        'MarkerSize',16,'LineWidth',1.5)
    else
        semilogy(xData,mean(Rend(indLow:indHigh,:),1),'.-',...
            'MarkerSize',16,'LineWidth',1.5,...
            'Color',colorList(iRange,:))
    end
    hold on
end
plot([0 xData(end)*2],[1 1],'k--','LineWidth',1.5)
if strcmp(paramToTest,'imageSizeCell')
    xlabel('Image width (px)')
else
    xlabel(paramToTest)
end
ylabel('R (%)')
legend(arrayfun(@(c) ['t = ' num2str(c,3) ' nm'],...
    tMean,'UniformOutput',false))

%% Plot % error in each peak for a certain setting

testToShow = 4;

indsToPlot = 1:7;
intDiff = abs(intValsArray(indsToPlot,:,testToShow+1)...
    -intValsArray(indsToPlot,:,testToShow));
pctE = 100*intDiff./intValsArray(indsToPlot,:,testToShow);

figure;
colorList = jet(nPlots).*0.8;
for iPlot = 1:nPlots
    semilogy(tArray,pctE(iPlot,:),...
        'Color',colorList(iPlot,:),'LineWidth',1.5)
    hold on
end
xlabel('Thickness (nm)')
ylabel('(I_{k+1} - I_{k}) / I_{k} (%)')
title(legendStr{testToShow})
legend(peakNames{indsToPlot})
