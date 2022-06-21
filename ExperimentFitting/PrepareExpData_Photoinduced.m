%% Prepare experimental diffraction peak change data

%% Load experimental data

fnameList = {'HWP18_FitGaussian_v2.mat',...
    'HWP20_FitGaussian_v2.mat',...
    'HWP22_FitGaussian_v2.mat',...
    'HWP25_FitGaussian_v2.mat',...
    'HWP27_FitGaussian_v2.mat',...
    'HWP29_FitGaussian_v2.mat',...
    'HWP31_FitGaussian_v2.mat',...
    'HWP33_FitGaussian_v2.mat',...
    'HWP35_FitGaussian_v2.mat'};
nFiles = length(fnameList);

% Timing parameters
tLim = 70; % Cut off at +70 ps, sample starts recovering after that 
% (or could be drift of the pump off the sample)
d0 = 25.6; % mm
psPermm = 6.67;

% HWP to fluence
HWPMeasArray = [5 10 15 20 25 30 35 40 45 50 55 60];
fluenceMeasArray = [1.15 0.124 0.257 1.57 3.91 6.99...
    10.4 13.7 16.8 18.9 19.9 19.9];
% fit a sin^2 wave
p0 = [max(fluenceMeasArray) 100 12]; %Amplitude, period, zero position
fitFunc = @(c,x) c(1).*sin(2*pi*((x-c(3))/c(2))).^2;
pFit = lsqcurvefit(fitFunc,p0,HWPMeasArray,fluenceMeasArray);
% figure;
% hold on
% plot(HWPMeasArray,fluenceMeasArray,'k.','MarkerSize',16);
% plot(0:60,fitFunc(cFit,0:60),'b-','LineWidth',1.5);
% xlabel('Half waveplate position')
% ylabel('Fluence (mJ/cm^2)')

% Compute fluence values for HWP used
HWPArray = cellfun(@(x) str2num(x(4:5)),fnameList);
fluenceArray = round(fitFunc(pFit,HWPArray),1);
plotLabels = cellfun(@(x) [num2str(x) ' mJ/cm^2'],...
    num2cell(fluenceArray),'UniformOutput',false);


setNameList = {'peakSet_200','peakSet_220','peakSet_400',...
    'peakSet_420','peakSet_440','peakSet_600','peakSet_620'};
setLabels = {'200','220','400','420','440','600','620'};
nSets = numel(setNameList);
hklExp = [2 0 0; 2 2 0; 4 0 0 ; 4 2 0; 4 4 0; 6 0 0; 6 2 0];

%% Extract mlogI for exp data

indsToAvgStart = 1:4;
indsToAvgEnd = 7:13;
nPts = 13;

IratioCell = cell(nFiles,1);
IfitCell = cell(nFiles,1);
Iratio = zeros(nSets,nPts);
mlogIStart = zeros(nFiles,nSets);
mlogIEnd = zeros(nFiles,nSets);
mlogIStartErr = zeros(nFiles,nSets);
mlogIEndErr = zeros(nFiles,nSets);
cErrArray = cell(nFiles,1);

for iFile = 1:nFiles
    load(fnameList{iFile});
    
    % Compute normalization skew, apply later to fit intensities
    normOffStart = mean(AnalysisDataOff.(setNameList{1}).normFactorArray(indsToAvgStart));
    normOnStart = mean(AnalysisDataOn.(setNameList{1}).normFactorArray(indsToAvgStart));
    normSkewStart = mean(normOnStart./normOffStart);
    
    normOffEnd = mean(AnalysisDataOff.(setNameList{1}).normFactorArray(indsToAvgEnd));
    normOnEnd = mean(AnalysisDataOn.(setNameList{1}).normFactorArray(indsToAvgEnd));
    normSkewEnd = mean(normOnEnd./normOffEnd);
    
    IratioStart = zeros(1,nSets);
    IratioEnd = zeros(1,nSets);
    for iSet = 1:nSets
        Ioff = mean(AnalysisDataOff.(setNameList{iSet}).volumeArrayNorm(1:nPts,:),2);
        IoffStart = mean(Ioff(indsToAvgStart));
        IoffEnd = mean(Ioff(indsToAvgEnd));
        Ion = mean(AnalysisDataOn.(setNameList{iSet}).volumeArrayNorm(1:nPts,:),2);
        IonStart = mean(Ion(indsToAvgStart));
        IonEnd = mean(Ion(indsToAvgEnd));
        Iratio(iSet,:) = Ion(1:nPts)./Ioff(1:nPts);
        IratioStart(iSet) = IonStart/IoffStart;
        IratioEnd(iSet) = IonEnd/IoffEnd;
    end
    IratioCell{iFile} = Iratio;
    
    dStagePos = scanArray{1}(1:nPts);
    
    tData = (dStagePos-d0)*psPermm;
    c0 = [IratioStart-IratioEnd, IratioStart];
    fitFun = @(c,x) multiExpConvGauss(...
        [0.5 0 c(1:nSets) 3.32*ones(1,nSets) c(nSets+1:end)],x,nSets);
    [cFit,~,residual,~,~,~,J] = lsqcurvefit(fitFun,c0,tData,Iratio(:));
    
    % Compute variance in residual for each plot
    residualMat = reshape(residual,[nSets,nPts]);
    varResPlots = var(residualMat,0,2);
    varResTotal = var(residual);
    
    preCov = inv(J'*J); % Without variance in residual added yet
    cErr = sqrt(diag(preCov).*...
        [varResPlots; varResPlots])';
    cErrArray{iFile} = cErr;
    
    figure;
    colorList = jet(nSets).*0.8;
    nFitPts = 200;
    tFit = linspace(tData(1),tData(end),nFitPts);
    Ifit = reshape(fitFun(cFit,tFit),[nSets,nFitPts]);
    hPlotList = gobjects(nSets,1);
    for iSet = 1:nSets
        hPlotList(iSet) = plot(tData,Iratio(iSet,:),'.',...
            'Color',colorList(iSet,:),'MarkerSize',12);
        hold on
        plot(tFit,Ifit(iSet,:),'-',...
            'Color',colorList(iSet,:));
    end
    legend(hPlotList,setLabels)
    xlabel('Time (ps)')
    ylabel('I_{on}/I_{off}')
    
    IfitCell{iFile} = Ifit;
    
    IratioStartFit = cFit(nSets+1:end).*normSkewStart;
    IratioStartErr = cErr(nSets+1:end).*normSkewStart;
    IratioEndFit = (cFit(nSets+1:end)-cFit(1:nSets)).*normSkewEnd;
    IratioEndErr = cErr(1:nSets).*normSkewEnd;
    
    mlogIStart(iFile,:) = -1*log(IratioStartFit);
    mlogIStartErr(iFile,:) = IratioStartErr./IratioStartFit;
    mlogIEnd(iFile,:) = -1*log(IratioEndFit);
    mlogIEndErr(iFile,:) = IratioEndErr./IratioEndFit;
    
end

colorList = jet(nFiles)*0.8;

% 420 peak plot vs fluence
figure;
hPlotList = gobjects(nFiles,1);
for iFile = 1:nFiles
    Iratio = IratioCell{iFile};
    hPlotList(iFile) = plot(tData,Iratio(4,:),'.',...
            'Color',colorList(iFile,:),'MarkerSize',10);
    hold on
    Ifit = IfitCell{iFile};
    plot(tFit,Ifit(4,:),'-',...
        'LineWidth',1.5,'Color',colorList(iFile,:));
end
legend(hPlotList,plotLabels)
xlabel('Time (ps)')
ylabel('I_{420}^{on}/I_{420}^{off}')
xlim([-20 70])

%% Save results

save('FluenceScan_mlogI.mat',...
    'fluenceArray',...
    'hklExp',...
    'mlogIStart','mlogIStartErr',...
    'mlogIEnd','mlogIEndErr')
