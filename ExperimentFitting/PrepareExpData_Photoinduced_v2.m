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
nStart = numel(indsToAvgStart);
indsToAvgEnd = 7:13;
nEnd = numel(indsToAvgEnd);
nPts = 13;

mlogIStart = zeros(nFiles,nSets);
mlogIEnd = zeros(nFiles,nSets);
mlogIStartErr = zeros(nFiles,nSets);
mlogIEndErr = zeros(nFiles,nSets);

for iFile = 1:nFiles
    load(fnameList{iFile});
    
    % Compute normalization skew, apply to Iratios and errors
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
        IoffStartErr = std(Ioff(indsToAvgStart))./sqrt(nStart);
        IoffEndErr = std(Ioff(indsToAvgEnd))./sqrt(nEnd);
        
        Ion = mean(AnalysisDataOn.(setNameList{iSet}).volumeArrayNorm(1:nPts,:),2);
        IonStart = mean(Ion(indsToAvgStart)).*normSkewStart;
        IonEnd = mean(Ion(indsToAvgEnd)).*normSkewEnd;
        IonStartErr = std(Ion(indsToAvgStart)).*normSkewStart./sqrt(nStart);
        IonEndErr = std(Ion(indsToAvgEnd)).*normSkewEnd./sqrt(nEnd);
        
        mlogIStart(iFile,iSet) = -log(IonStart/IoffStart);
        mlogIEnd(iFile,iSet) = -log(IonEnd/IoffEnd);
        mlogIStartErr(iFile,iSet) = IonStartErr./IonStart ...
            + IoffStartErr./IoffStart;
        mlogIEndErr(iFile,iSet) = IonEndErr./IonEnd ...
            + IoffEndErr./IoffEnd;
    end
    
    
%     dStagePos = scanArray{1}(1:nPts); 
%     tData = (dStagePos-d0)*psPermm;
    
end

colorList = jet(nSets)*0.8;

% Plot on, off mlogIratios and errors
figure;
hPlotList = gobjects(nSets,1);
for iSet = 1:nSets
    hPlotList(iSet) = errorbar(fluenceArray,...
        mlogIStart(:,iSet),mlogIStartErr(:,iSet),...
        '.','Color',colorList(iSet,:),'MarkerSize',10);
    hold on
end
legend(setLabels)
xlabel('Fluence (mJ cm^{-2})')
ylabel('-log(I_{on}/I_{off}) for t < -4 ps')

figure;
hPlotList = gobjects(nSets,1);
for iSet = 1:nSets
    hPlotList(iSet) = errorbar(fluenceArray,...
        mlogIEnd(:,iSet),mlogIEndErr(:,iSet),...
        '.','Color',colorList(iSet,:),'MarkerSize',10);
    hold on
end
legend(setLabels)
xlabel('Fluence (mJ cm^{-2})')
ylabel('-log(I_{on}/I_{off}) for t > 16 ps')

%% Save results

save('FluenceScan_mlogI_v2.mat',...
    'fluenceArray',...
    'hklExp',...
    'mlogIStart','mlogIStartErr',...
    'mlogIEnd','mlogIEndErr')
