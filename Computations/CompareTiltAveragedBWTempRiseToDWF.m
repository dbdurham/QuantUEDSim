%% Compare tilt-averaged BW at two temperatures to DW model
% Figure 2b-c in D.B. Durham et al, Arxiv preprint 2022

%% Load stack of BW simulated patterns at T1

[filename, pathname] = uigetfile('*.mat','Load DP library');
load([pathname filename])

sDiff1 = sDiff;
Ilib1 = Ilib;
iEnd = size(Ilib1,5);

%% Load stack of BW simulated patterns at T2

[filename, pathname] = uigetfile('*.mat','Load DP library');
load([pathname filename])

sDiff2 = sDiff;
Ilib2 = Ilib;
iEnd2 = size(Ilib2,5);

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

[GhklTest,GmagTest] = computeScatteringVectors(hklTest,sDiff.Gvec);

%% Plot the intensity change for a set tilt spread against s^2

iTheta = 40;

iUCToPlot = 1:10:100;
nPlots = numel(iUCToPlot);

IArraySim2 = extractIntsFromDP(Ilib2(:,:,:,iTheta,iEnd2),...
        sDiff2.qxaStore,sDiff2.qyaStore,GhklTest);
IArraySim1 = extractIntsFromDP(Ilib1(:,:,:,iTheta,iEnd),...
        sDiff1.qxaStore,sDiff1.qyaStore,GhklTest);
mlogIRatioArray = -log(IArraySim2./IArraySim1);

colorList = parula(nPlots);

sList = 2*pi*GmagTest;
s2 = 4*pi^2*GmagTest.^2;
du2BW = sDiff2.uRMS^2-sDiff1.uRMS^2;

% colorList = 0.8*jet(nPeaks);
% figure('Position',[200 200 400 300]);
% lobjs = gobjects(nPeaks,1);
% plot([0 tArray(end)],...
%         0.*[1 1],...
%         'k-','LineWidth',0.5)
% hold on
% for iPlot = 1:nPeaks
%     du2DW = mlogIRatioArray(iPlot,:)./sList(iPlot).^2;
%     pctErrorDWSim = 100*(du2DW-du2BW)./du2BW;
%     lobjs(iPlot) = plot(tArray,...
%         pctErrorDWSim,...
%         '-','LineWidth',1.5,...
%         'Color',colorList(iPlot,:));
%     
% end
% xlim([0 40])
% xlabel('Thickness(nm)')
% ylabel('\DeltaT_{Kin-Dyn} (%)')
% legend(lobjs,peakNames)

colorList = 0.8*jet(nPeaks);
figure;
lobjs = gobjects(nPeaks,1);
pctChangeArray = 100*(IArraySim2-IArraySim1)./IArraySim1;
pctChangeDW = 100*(exp(-s2*du2BW)-1);
for iPlot = 1:nPeaks
    plot([0 tArray(end)],...
        pctChangeDW(iPlot).*[1 1],...
        '--','LineWidth',0.75,...
        'Color',colorList(iPlot,:))
    hold on
    lobjs(iPlot) = plot(tArray,...
        pctChangeArray(iPlot,:),...
        '-','LineWidth',1.5,...
        'Color',colorList(iPlot,:));
    
end
xlabel('Thickness (nm)')
ylabel('Intensity change (%)')
legend(lobjs,peakNames)
title(['\Delta u^2 = ' num2str(3*du2BW,3) ...
    ', \sigma_{\theta} = ' num2str(sigmaThetaSamp(iTheta)*1e3) ' mrad'])

%% Generate DW fit map 

sList = 2*pi*GmagTest;
du2BW = sDiff2.uRMS^2-sDiff1.uRMS^2;

nUC = size(Ilib2,3);

du2Array = zeros(nUC,nTheta);
for iTheta = 1:nTheta
    IArraySim2 = extractIntsFromDP(Ilib2(:,:,:,iTheta,iEnd2),...
        sDiff2.qxaStore,sDiff2.qyaStore,GhklTest);
    IArraySim1 = extractIntsFromDP(Ilib1(:,:,:,iTheta,iEnd),...
        sDiff1.qxaStore,sDiff1.qyaStore,GhklTest);
    mlogIRatioArray = -log(IArraySim2./IArraySim1);
    du2Array(:,iTheta) = sum(mlogIRatioArray.*sList.^2)...
                ./sum(sList.^4);
end

du2ErrorArray = 100*(du2Array-du2BW)./du2BW;

% cpts = [0   0   0.35 0;
%         0   0   1   0.36;
%         1   1   1   0.5;
%         1   0   0   0.63;
%         0.2 0   0   1];
% cmap = generateGradColormap(cpts,nColors);

nColors = 1024;
cpts = [0.5   1   0.8   0;
        0.3   0.4   0.8   0.36;
        0   0   0   0.5;
        1   0   0   0.63;
        1   1   0   1];
cmap = generateGradColormap(cpts,nColors);

contourColor1 = 0.8.*[204 164 0]./255;
contourColor2 = 0.85.*[204 164 0]./255;
contourColor3 = [204 164 0]./255;

figure;
imagesc(sigmaThetaSamp([1 end])*1e3,...
    tArray([1 end]),...
    du2ErrorArray);
hold on
[X,Y] = meshgrid(sigmaThetaSamp.*1e3,tArray);
% contour(X,Y,du2ErrorArray,[-10,-10],'--','LineWidth',1,'Color',contourColor1)
contour(X,Y,du2ErrorArray,[-25,-25],'--','LineWidth',1,'Color','w')
contour(X,Y,du2ErrorArray,[-50,-50],'--','LineWidth',1,'Color','w')
contour(X,Y,du2ErrorArray,[-75,-75],'--','LineWidth',1,'Color','w')
xlabel('\sigma_{\theta} (mrad)')
ylabel('Thickness (nm)')
title('Error in \Deltau^2 using DW fit (%)')
colormap(cmap)
colorbar()
caxis([-250 250])
set(gca,'ydir','normal')
xlim([1.25 125])

