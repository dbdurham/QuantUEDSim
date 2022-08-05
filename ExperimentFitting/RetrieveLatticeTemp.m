%% Fit bloch wave results to extract lattice temperature

%% Define and load experiment data

load('FluenceScan_mlogI_v2.mat')

nPeaks = size(hklExp,1);

peakNames = cell(nPeaks,1);
for iPeak = 1:nPeaks
    peakNames{iPeak} = strrep(num2str(hklExp(iPeak,:)),' ','');
end

nFilesExp = numel(fluenceArray);

plotLabels = cellfun(@(x) [num2str(x) ' mJ/cm^2'],...
    num2cell(fluenceArray),'UniformOutput',false);

%% Load and prepare bloch wave library for matching

% Fix thickness, angle spread to best fit for static pattern
iUC = 33; 
iTheta = 38; 

% Define simulated library to use for matching
ftag = ['C:/Users/ddurh/OneDrive/Documents/Berkeley/Research'...
    '/HiRES/DiffractionSimPackageOutputs/TestOutputs/Library_u2Tilt_BW_01/'...
    '750keV_160mrad_50UC_BW_u'];	
fnameList = arrayfun(@(c) [ftag num2str(c,'%02.f') '.mat'],...
    1:15,'UniformOutput',false);

load(fnameList{1})
nFiles = numel(fnameList);

iEnd = size(Ilib,5);
tArray = (1:nUC)*sDiff.cellDim(3)*0.1; %nm

[GhklExp,GmagExp] = computeScatteringVectors(hklExp,sDiff.Gvec);
sList = 2*pi*GmagExp';

uRMS = zeros(nFiles,1);
intsSim = zeros(nPeaks,nFiles);

for iFile = 1:nFiles
    load(fnameList{iFile})
    
    % Extract peak intensities vs thickness   
    intsSim(:,iFile) = extractIntsFromDP(Ilib(:,:,iUC,iTheta,iEnd),...
        sDiff.qxaStore,sDiff.qyaStore,GhklExp);
    
    % Store uRMS
    uRMS(iFile) = sDiff.uRMS;

end

mlogIArraySim = -log(intsSim./intsSim(:,1));
du2Sim = uRMS.^2 - uRMS(1).^2; % 1D mean-square displacements

% Interpolate w.r.s. du2
mlogISim = @(du2) interp1(du2Sim,mlogIArraySim',du2,'pchip');

% Visualize model library
du2Plot = linspace(0,du2Sim(end),100);
mlogIPlot = mlogISim(du2Plot);
figure;
colorList = 0.8.*jet(nPeaks);
lobjs = gobjects(nPeaks,1);
for iPeak = 1:nPeaks
    lobjs(iPeak) = plot(3*du2Plot,mlogIPlot(:,iPeak),'-',...
        'LineWidth',1.5,...
        'Color',colorList(iPeak,:));
    hold on
    plot(3*du2Sim,mlogIArraySim(iPeak,:),'.',...
        'MarkerSize',12,...
        'Color',colorList(iPeak,:));
end
legend(lobjs,peakNames)
xlabel('\Deltau^2')
ylabel('-log(I_{on}/I_{off})')

%% Perform DW and MS fits

DWFitType = 'LSQ'; % Choices are LSQ, Chi2

du2Start = zeros(nFilesExp,1);
du2End = zeros(nFilesExp,1);
du2SimStart = zeros(nFilesExp,1);
du2SimEnd = zeros(nFilesExp,1);

for iFile = 1:nFilesExp
    
    % DW fit
    switch DWFitType
        case 'Chi2'            
            [du2Start(iFile),~] = linRegMLENoInt(...
                sList.^2,...
                mlogIStart(iFile,:),...
                mlogIStartErr(iFile,:));
            [du2End(iFile),~] = linRegMLENoInt(...
                sList.^2,...
                mlogIEnd(iFile,:),...
                mlogIEndErr(iFile,:));
        case 'LSQ'
            du2Start(iFile) = sum(mlogIStart(iFile,:).*sList.^2)...
                ./sum(sList.^4);
            du2End(iFile) = sum(mlogIEnd(iFile,:).*sList.^2)...
                ./sum(sList.^4);
    end
    
    % Simulation fit
    residualSim = @(p,mlogIExp) sum(abs(mlogISim(p) - mlogIExp).^2);
    du2SimStart(iFile) = fminunc(...
        @(c) residualSim(c,mlogIStart(iFile,:)),du2Start(iFile));
    du2SimEnd(iFile) = fminunc(...
        @(c) residualSim(c,mlogIEnd(iFile,:)),du2End(iFile));
end

%% Visualize fit results

colorList = [0 0 0;...
    154 22 163; ...
    177 26 88; ...
    208 30 38; ...
    209 98 33; ...
    212 150 28; ...
    216 207 27; ...
    218 213 123; ...
    220 220 220]./255;
    

% Show fits before time zero
figure;
hPlotList = gobjects(nFilesExp,1);
for iFile = 1:nFilesExp
    switch DWFitType
        case 'Chi2'
            hPlotList(iFile) = errorbar(...
                sList.^2,mlogIStart(iFile,:),mlogIStartErr(iFile,:),...
                '.','MarkerSize',14,'Color',colorList(iFile,:));
        case 'LSQ'
            hPlotList(iFile) = plot(...
                sList.^2,mlogIStart(iFile,:),...
                '.','MarkerSize',14,'Color',colorList(iFile,:));
    end
    hold on
    plot([0 100],[0 100*du2Start(iFile)],'-',...
        'MarkerSize',14,'Color',colorList(iFile,:));
    plot(sList.^2,mlogISim(du2SimStart(iFile)),...
        'd','MarkerSize',10,'Color',colorList(iFile,:))
end
legend(hPlotList,plotLabels)
xlabel('s^{2} (Angstroms^{-1})')
ylabel('-ln(I_{on}/I_{off}) before t_0')
xlim([0 100])

% Show fits after time zero
figure('Position',[200 200 400 300]);
hPlotList = gobjects(nFilesExp,1);
for iFile = 1:nFilesExp
    switch DWFitType
        case 'Chi2'
            hPlotList(iFile) = errorbar(...
                sList.^2,mlogIEnd(iFile,:),mlogIEndErr(iFile,:),'.',...
                'MarkerSize',14,'Color',colorList(iFile,:));
        case 'LSQ'
            hPlotList(iFile) = plot(...
                sList.^2,mlogIEnd(iFile,:),...
                '.','MarkerSize',14,'Color',colorList(iFile,:));
    end
    hold on
    plot([0 100],[0 100*du2End(iFile)],'-',...
        'MarkerSize',14,'Color',colorList(iFile,:));
    plot(sList.^2,mlogISim(du2SimEnd(iFile)),...
        'd','MarkerSize',10,'Color',colorList(iFile,:))
end
legend(hPlotList,plotLabels)
xlabel('s^{2} (Angstroms^{-1})')
ylabel('-ln(I_{on}/I_{off}) after t_0')
xlim([0 100])
% Show one after time zero fit
iFile = 6;
q2 = (sList.^2)/(4*pi^2);
figure('Position',[200 200 400 300]);
hPlotList = gobjects(nFilesExp,1);
errorbar(...
    q2,mlogIEnd(iFile,:),mlogIEndErr(iFile,:),'k.',...
    'MarkerSize',14);
hold on
plot([0 100/(4*pi^2)],[0 100*du2End(iFile)],'-',...
    'MarkerSize',14,'Color',[0.8 0.6 0]);
plot(q2,mlogISim(du2SimEnd(iFile)),...
    'd','MarkerSize',8,'Color',[0 0.6 0])
legend('Exp','Kin fit','Dyn fit')
xlabel('q^{2} (Angstroms^{-1})')
ylabel('-ln(I_{on}/I_{off}) after t_0')
xlim([0 100/(4*pi^2)])


% Scatter plot comparing agreement for max fluence after time zero
figure;
plotLims = [-0.05 0.3];
plot(mlogIEnd(end,:),du2End(end).*sList.^2,'r.','MarkerSize',12)
hold on
plot(mlogIEnd(end,:),mlogISim(du2SimEnd(end)),'b.','MarkerSize',12)
plot(plotLims,plotLims,'k--','LineWidth',1.5)
xlabel('-ln(I_{on,exp}/I_{off,exp})')
ylabel('-ln(I_{on,sim}/I_{off,sim})')
legend('DW','MS')

% Compute and display residuals
resValSim = zeros(1,nFiles);
resValDW = zeros(1,nFiles);
for iFile = 1:nFiles
    resValSim(iFile) = sum(abs(mlogISim(du2SimEnd(iFile)) - mlogIEnd(iFile,:)).^2)...
        ./sum(abs(mlogIEnd(iFile,:)).^2);
    resValDW(iFile) = sum(abs(du2End(iFile).*sList.^2 - mlogIEnd(iFile,:)).^2)...
        ./sum(abs(mlogIEnd(iFile,:)).^2);
end
figure;
plot(fluenceArray,100*resValSim,'kd',...
    'MarkerSize',6,'MarkerFaceColor','k')
hold on
plot(fluenceArray,100*resValDW,'k.','MarkerSize',18)
xlabel('Peak fluence (mJ/cm^{2})')
ylabel('Residual (%)')

% Thermal motions vs fluence (3D)
figure;
plot(fluenceArray,3*du2Start,'k.','MarkerSize',18)
hold on
plot(fluenceArray,3*du2SimStart,'kd',...
    'MarkerSize',6,'MarkerFaceColor','k')
plot(fluenceArray,3*du2End,'r.','MarkerSize',18)
plot(fluenceArray,3*du2SimEnd,'rd',...
    'MarkerSize',6,'MarkerFaceColor','r')
xlabel('Peak fluence (mJ/cm^{2})')
ylabel('\Delta u_{rms,3D}^{2} (Angstroms^{2})')
legend('DW before t_0','Sim before t_0',...
    'DW after t_0','Sim after t_0')
xlim([0 11])

%% Compute and add lattice temperature to curves

% Lattice temperature vs fluence (using model proposed by 
% Owens and Williams)
h = 6.626e-34; % Planck's constant, J*s
m = 196.97*(1.662e-27); % Gold atomic mass, kg
k = 1.38e-23; % Boltzmann constant, J/K
TDebye = 175; % Debye Temp for gold at 300K, K

gam = 3.03; % Gruneisen constant
alpha = 13.9e-6; % Linear coef of thermal expansion for gold, 1/K

T0 = 300;

Y = @(T) -((12*h^2)/(m*k*TDebye^2)) * (T-T0) ...
    .* (1 + 2*alpha*gam*T) *1e20; % Angstroms^2

du2 = @(T) -Y(T)/(16*pi^2); % Angstroms^2

% Linearize over 300 K range:
Kperu2 = 300/(du2(600)-du2(300));

figure;
plotLims = [-0.001 0.014];
plot(fluenceArray,3*du2Start,'k.','MarkerSize',18)
hold on
plot(fluenceArray,3*du2SimStart,'kd',...
    'MarkerSize',6,'MarkerFaceColor','k')
plot(fluenceArray,3*du2End,'r.','MarkerSize',18)
plot(fluenceArray,3*du2SimEnd,'rd',...
    'MarkerSize',6,'MarkerFaceColor','r')
xlabel('Peak fluence (mJ/cm^{2})')
ylabel('\Delta u_{rms,3D}^{2} (Angstroms^{2})')
legend('DW before t_0','MS before t_0',...
    'DW after t_0','MS after t_0')
xlim([0 11])
ylim(plotLims)
yyaxis right
ylim(plotLims./3*Kperu2)
ylabel('\Delta T (K)')
set(gca,'YColor',[0 0 0])


figure('Position',[200 200 380 300]);
plotLims = [0 0.011];
xLims = [0 11];
dTPerF = 11.61; % Predicted
plot(xLims,xLims*dTPerF/(Kperu2/3),'k--','LineWidth',1.5);
hold on
plot(fluenceArray,3*(du2End-du2Start),'r.','MarkerSize',18)
hold on
plot(fluenceArray,3*(du2SimEnd-du2SimStart),'bd',...
    'MarkerSize',6,'MarkerFaceColor','b')
xlabel('Peak Fluence (mJ cm^{-2})')
ylabel('\Delta u^{2} (Angstroms^{2})')
xlim(xLims)
ylim(plotLims)
yyaxis right
ylim(plotLims./3*Kperu2)
ylabel('\Delta T (K)')
set(gca,'YColor',[0 0 0])

% fit line to data
mFit = sum(3*(du2End-du2Start).*fluenceArray')...
    ./sum(fluenceArray.^2);
mFitSim = sum(3*(du2SimEnd-du2SimStart).*fluenceArray')...
    ./sum(fluenceArray.^2);

%plot(xLims,xLims*mFit,'k-','LineWidth',1.5);
legend('Predicted','DW','Sim')