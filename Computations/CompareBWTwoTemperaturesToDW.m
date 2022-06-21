%% Compare thickness dependence: BW vs MS

% Set up Bloch Wave simulations
options1.uRMS = 0.0894;
options2.uRMS = 0.1033;
sDiff1 = setupSimBW(options1);
sDiff2 = setupSimBW(options2);

% Inputs for calculation
theta1 = 0; % rad, x component of tilt
theta2 = 0; % rad, y component of tilt
nUC = 100; % Number of unit cells to simulate

IDiff1 = calcDiffBW(-theta1,-theta2,nUC,sDiff1);
IDiff2 = calcDiffBW(-theta1,-theta2,nUC,sDiff2);

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

[GhklTest,GmagTest] = computeScatteringVectors(hklTest,sDiff1.Gvec);

IArray1 = extractIntsFromDP(IDiff1,...
    sDiff1.qxaStore,sDiff1.qyaStore,GhklTest);
IArray2 = extractIntsFromDP(IDiff2,...
    sDiff2.qxaStore,sDiff2.qyaStore,GhklTest);
pctChangeArray = 100*(IArray2-IArray1)./IArray1;
mlogIRatioArray = -log(IArray2./IArray1);

tArray = (1:nUC)*sDiff1.cellDim(3)*0.1; %nm

sList = 2*pi*GmagTest;
du2BW = sDiff2.uRMS^2-sDiff1.uRMS^2;

colorList = 0.8*jet(nPeaks);
figure('Position',[200 200 400 300]);
lobjs = gobjects(nPeaks,1);
plot([0 tArray(end)],...
        0.*[1 1],...
        'k-','LineWidth',0.5)
hold on
for iPlot = 1:nPeaks
    du2DW = mlogIRatioArray(iPlot,:)./sList(iPlot).^2;
    pctErrorDWSim = 100*(du2DW-du2BW)./du2BW;
    lobjs(iPlot) = plot(tArray,...
        pctErrorDWSim,...
        '-','LineWidth',1.5,...
        'Color',colorList(iPlot,:));
    
end
xlabel('Thickness(nm)')
ylabel('\DeltaT_{Kin-Dyn} (%)')
legend(lobjs,peakNames)

% colorList = 0.8*jet(nPeaks);
% figure('Position',[200 200 400 300]);
% lobjs = gobjects(nPeaks,1);
% plot([0 tArray(end)],...
%         0.*[1 1],...
%         'k-','LineWidth',0.5)
% hold on
% plot([0 tArray(end)],...
%         -du2BW.*[1 1],...
%         'k--','LineWidth',0.5)
% for iPlot = 1:nPeaks
%     du2DW = mlogIRatioArray(iPlot,:)./sList(iPlot).^2;
%     lobjs(iPlot) = plot(tArray,...
%         -du2DW,...
%         '-','LineWidth',1.5,...
%         'Color',colorList(iPlot,:));
%     
% end
% xlabel('Thickness(nm)')
% ylabel('log(I_{on}/I_{off})/s^2')
% legend(lobjs,peakNames)
% xlim([0 40])

colorList = 0.8*jet(nPeaks);
figure('Position',[200 200 400 300]);
lobjs = gobjects(nPeaks,1);
plot([0 tArray(end)],...
        0.*[1 1],...
        'k-','LineWidth',0.5)
hold on
for iPlot = 1:nPeaks  
    plot([0 tArray(end)],...
            -du2BW.*s2(iPlot).*[1 1],...
            '--','LineWidth',0.5,...
            'Color',colorList(iPlot,:))
    du2DW = mlogIRatioArray(iPlot,:);
    lobjs(iPlot) = plot(tArray,...
        -du2DW,...
        '-','LineWidth',1.5,...
        'Color',colorList(iPlot,:));
    
end
xlabel('Thickness(nm)')
ylabel('I_{on}/I_{off}')
legend(lobjs,peakNames)
xlim([0 40])

colorList = 0.8*jet(nPeaks);
figure;
lobjs = gobjects(nPeaks,1);
pctChangeDW = 100*(exp(-sList.^2*du2BW)-1);
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
title(['\Delta u^2 = ' num2str(3*du2BW,3)])

colorList = 0.8*jet(nPeaks);
figure;
lobjs = gobjects(nPeaks,1);
pctChangeDW = 100*(exp(-sList.^2*du2BW)-1);
for iPlot = 1:nPeaks
    plot([0 tArray(end)],...
        sign(pctChangeDW(iPlot)).*log10(pctChangeDW(iPlot)).*[1 1],...
        '--','LineWidth',0.75,...
        'Color',colorList(iPlot,:))
    hold on
    lobjs(iPlot) = plot(tArray,...
        sign(pctChangeArray(iPlot,:)).*log10(pctChangeArray(iPlot,:)),...
        '-','LineWidth',1.5,...
        'Color',colorList(iPlot,:));
    
end
xlabel('Thickness (nm)')
ylabel('Intensity change (%)')
legend(lobjs,peakNames)
title(['\Delta u^2 = ' num2str(3*du2BW,3)])