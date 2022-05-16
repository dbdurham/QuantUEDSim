%% Compute a tilt series and compare the results between 
% parabolic and spherical propagators

%% 1. Set up simulation cell + parameters

% Can set parameters here through the options struct, 
% or in the setup function itself
options.E0 = 750e3; % Electron energy (keV)
options.uRMS = 0.0894; % 1D rms atomic displacement
options.qRange = 4; % q range to store in *output* (A^-1)
options.cellMult = 1;
sDiff = setupMultisliceSim(options); 

%% 2. Compute tilt series of diffraction patterns

sDiff.storeRealSpaceCellEWs = false;
nUCs = 100; % Number of simulation cells to propagate through

thetaList = [0:0.0025:0.05 0.06:0.01:0.4];
nTheta = numel(thetaList);

for iTheta = 1:nTheta
    disp(['Computing tilt #' num2str(iTheta)])
    theta_x = thetaList(iTheta); % Incident angle of the e-beam (rad)
    theta_y = 0.0;
    coefs = [nUCs,0,theta_x,theta_y];
    [~,EWDiff,~,~] = calcDiff(sDiff,coefs);
    IDiff = double(abs(EWDiff).^2);
    
    if iTheta == 1
        IDiffStack = zeros([size(IDiff) nTheta]);
    end
    IDiffStack(:,:,:,iTheta) = IDiff;
end

%% 2b. Visualize diffraction peaks

% list q values of peaks to study
h1 = [0.4902 0]; % basis vectors
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

% Determine peak locations   
indPeaks = zeros(nPeaks,1);
for iPeak = 1:nPeaks
    [~,indPeaks(iPeak)] = ...
        min((qxyPeaks(iPeak,2)-sDiff.qxaStore(:)).^2 ...
        + (qxyPeaks(iPeak,1)-sDiff.qyaStore(:)).^2);
end

% Extract intensity values
intValsArray = zeros(nPeaks,nUCs,nTheta);
I0Array = zeros(nUCs,nTheta);
for iTheta = 1:nTheta
    for iUC = 1:nUCs
        I = IDiffStack(:,:,iUC,iTheta);
        I0 = I(1,1);
        I0Array(iUC,iTheta) = I0;
        % normalize peaks by total *diffracted* intensity
        % I = I./(sum(I(:))-I0); 
        intValsArray(:,iUC,iTheta) = I(indPeaks);
    end
end

tArray = 0.1*sDiff.cellDim(3)*(1:nUCs);

% Thickness dependence of a selected angle
iTheta = 1;

figure;
subplot(2,1,1)
plot(tArray,I0Array(:,iTheta),'k-','LineWidth',1.5)
xlabel('Thickness (nm)')
ylabel('Primary beam intensity')
xlim([0 tArray(end)])

subplot(2,1,2)
nPlots = nPeaks;
colorList = jet(nPlots).*0.8;
for iPlot = 1:nPlots
    plot(tArray,squeeze(intValsArray(iPlot,:,iTheta)),...
            '-','Color',colorList(iPlot,:),'LineWidth',1.5)
    hold on
end
xlabel('Thickness (nm)')
ylabel('Diffracted intensity')
legend(peakNames)
xlim([0 tArray(end)])

% Tilt angle dependence of a selected thickness
iUC = 30;

figure;
subplot(2,1,1)
plot(thetaList*1e3,I0Array(iUC,:),'k-','LineWidth',1.5)
xlabel('Sample tilt (mrad)')
ylabel('Primary beam intensity')
xlim([0 thetaList(end)*1e3])

subplot(2,1,2)
nPlots = nPeaks;
colorList = jet(nPlots).*0.8;
for iPlot = 1:nPlots
    plot(thetaList*1e3,squeeze(intValsArray(iPlot,iUC,:)),...
            '-','Color',colorList(iPlot,:),'LineWidth',1.5)
    hold on
end
xlabel('Sample tilt (mrad)')
ylabel('Diffracted intensity')
legend(peakNames)
xlim([0 thetaList(end)*1e3])

% Save analysis

[savename,savepath] = uiputfile('*.mat');
save([savepath,savename],'sDiff','intValsArray')