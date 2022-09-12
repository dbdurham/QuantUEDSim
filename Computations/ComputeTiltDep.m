%% Compute tilt dependence of diffracted intensities

%% 1. Set up simulation cell + parameters

% Can set parameters here through the options struct, 
% or in the setup function itself
simType = 'Multislice'; %Simulation method: 
% Kinematical, Bloch Waves, or Multislice
options = struct;
options.E0 = 750e3; % Electron energy (eV)
options.uRMS = 0.0894; % 1D rms atomic displacement
options.qRange = 4; % q range to store in *output* (A^-1)
options.cellMult = 1;
options.downSampFacCell = 1;
sDiff = setupSim(simType,options); 

%% 2. Compute tilt series of diffraction patterns

sDiff.storeRealSpaceCellEWs = false;
nUCs = 100; % Number of simulation cells to propagate through

thetaList = [0:0.0025:0.05 0.06:0.01:0.4]; %rad
nTheta = numel(thetaList);

for iTheta = 1:nTheta
    disp(['Computing tilt #' num2str(iTheta)])
    theta1 = thetaList(iTheta); % Incident angle of the e-beam (rad)
    theta2 = 0.0;
    IDiff = calcDiff(theta1,theta2,nUCs,sDiff);
    
    if iTheta == 1
        IDiffStack = zeros([size(IDiff) nTheta]);
    end
    IDiffStack(:,:,:,iTheta) = IDiff;
end

%% 2b. Visualize diffraction peaks

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

% Extract intensity values
IArray = zeros(nPeaks,nUCs,nTheta);
I0Array = zeros(nUCs,nTheta);
for iTheta = 1:nTheta
    IDiff = IDiffStack(:,:,:,iTheta);
    IArray(:,:,iTheta) = extractIntsFromDP(IDiff,...
        sDiff.qxaStore,sDiff.qyaStore,GhklTest);
    I0Array(:,iTheta) = extractIntsFromDP(IDiff,...
        sDiff.qxaStore,sDiff.qyaStore,[0 0 0]);
end

tArray = 0.1*sDiff.cellDim(3)*(1:nUCs);

% Thickness dependence of a selected angle
iTheta = 1;

showIvt(IArray(:,:,iTheta),I0Array(:,iTheta),tArray,peakNames);
title(['\theta_1 = ' num2str(1e3*thetaList(iTheta)) ' mrad'])

% Tilt angle dependence of a selected thickness
iUC = 30;

figure;
subplot(2,1,1)
plot(thetaList*1e3,I0Array(iUC,:),'k-','LineWidth',1.5)
xlabel('Sample tilt (mrad)')
ylabel('Primary beam intensity')
xlim([0 thetaList(end)*1e3])
title(['Thickness = ' num2str(tArray(iUC)) ' nm'])

subplot(2,1,2)
nPlots = nPeaks;
colorList = jet(nPlots).*0.8;
for iPlot = 1:nPlots
    plot(thetaList*1e3,squeeze(IArray(iPlot,iUC,:)),...
            '-','Color',colorList(iPlot,:),'LineWidth',1.5)
    hold on
end
xlabel('Sample tilt (mrad)')
ylabel('Diffracted intensity')
legend(peakNames)
xlim([0 thetaList(end)*1e3])

%% Save analysis

[savename,savepath] = uiputfile('*.mat');
save([savepath,savename],'sDiff','IArray','I0Array')