%% Compare temperature effect on peak intensities to Debye-Waller Factor
%% in a flat crystal
% as used in D.B. Durham et al, ArXiv Preprint 2022 -- Fig 4a

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
lobjs = gobjects(nPeaks,1);
pctChangeDW = 100*(exp(-sList.^2*du2BW)-1);
axesLimits = {[0 4],[-2 0]};
yTicks = {[0 1 2 3 4],[-2 -1 0]};
yTickLabels = {{'10^0','10^1','10^2','10^3','10^4'},...
    {'-10^{2}','-10^{1}','-10^0'}};
subPos = {[],[]};
inRange = {pctChangeArray>0,pctChangeArray<0};
for iSub = 1:2
    subplot(2,1,iSub)
    for iPlot = 1:nPeaks
        yData = log10(pctChangeArray(iPlot,:));
        yData(~inRange{iSub}(iPlot,:)) = -1;
        if iSub == 2
            yData = yData.*-1;
        end
        lobjs(iPlot) = plot(tArray,...
            yData,...
            '-','LineWidth',1.5,...
            'Color',colorList(iPlot,:));
        hold on
        if iSub == 2
            plot([0 tArray(end)],...
                -1*log10(pctChangeDW(iPlot)).*[1 1],...
                '--','LineWidth',0.75,...
                'Color',colorList(iPlot,:))
        end
    end
    if iSub == 1
        xticks([])
        legend(lobjs,peakNames)
        title(['\Delta u^2 = ' num2str(3*du2BW,3)])
    end
    xlim([0 40])
    ylim(axesLimits{iSub})
    yticks(yTicks{iSub})
    yticklabels(yTickLabels{iSub})
end
xlabel('Thickness (nm)')
ylabel('Intensity change (%)')
