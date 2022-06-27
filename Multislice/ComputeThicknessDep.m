%% Compute thickness dependent diffraction peaks + wave functions

%% 1. Set up simulation cell + parameters

% Can set parameters here through the options struct, 
% or in the setup function itself
options.E0 = 750e3; % Electron energy (eV)
options.uRMS = 0.0894; % 1D rms atomic displacement
options.qRange = 4; % q range to store in *output* (A^-1)
options.cellMult = 1;
sDiff = setupMultisliceSim(options); 

%% 2a. Compute diffraction patterns
sDiff.storeRealSpaceCellEWs = false;
nUCs = 100; % Number of simulation cells to propagate through
theta_x = 0.0; % Incident angle of the e-beam (rad)
theta_y = 0.0;
coefs = [nUCs,0,theta_x,theta_y];
[~,EWDiff,~] = calcDiff(sDiff,coefs);
IDiff = double(abs(EWDiff).^2);

tArray = 0.1*sDiff.cellDim(3)*(1:nUCs);
StackViewerDiff(fftshift(fftshift(IDiff,1),2),tArray)

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
intValsArray = zeros(nPeaks,nUCs);
I0Array = zeros(1,nUCs);
for iUC = 1:nUCs
    I = IDiff(:,:,iUC);
    I0 = I(1,1);
    I0Array(iUC) = I0;
    % normalize peaks by total *diffracted* intensity
    % I = I./(sum(I(:))-I0); 
    intValsArray(:,iUC) = I(indPeaks);
end

tArray = 0.1*sDiff.cellDim(3)*(1:nUCs); % Thickness (nm)

figure;
subplot(2,1,1)
plot(tArray,I0Array,'k-','LineWidth',1.5)
xlabel('Thickness (nm)')
ylabel('Primary beam intensity')
xlim([0 tArray(end)])

subplot(2,1,2)
nPlots = nPeaks;
colorList = jet(nPlots).*0.8;
for iPlot = 1:nPlots
    plot(tArray,squeeze(intValsArray(iPlot,:)),...
            '-','Color',colorList(iPlot,:),'LineWidth',1.5)
    hold on
end
xlabel('Thickness (nm)')
ylabel('Diffracted intensity')
legend(peakNames)
xlim([0 tArray(end)])

%% 2c. Extract kinematical diff signals

fname = ['../2022-01-16_RippledKinematicalSims/' ... 
    '2022-01-29_KinSims0Deg/KinDiff_Gold_60keV_u0p0894_0Deg.mat'];
load(fname,'Ghkl','I');

nPeaks = size(qxyPeaks,1);
indPeaks = zeros(nPeaks,1);
intKinArray = zeros(nPeaks,nUCs);
for iPeak = 1:nPeaks
    [~,indPeaks(iPeak)] = ...
        min((qxyPeaks(iPeak,1)-Ghkl(:,1)).^2 ...
        + (qxyPeaks(iPeak,2)-Ghkl(:,2)).^2);
    intKinArray(iPeak,:) = I(indPeaks(iPeak),:);
end



%% 2d. Overlay kinematical and multislice signals

figure('Position',[100 100 400 250]);
nPlots = nPeaks;
colorList = jet(nPlots).*0.8;
plotObjs = gobjects(nPlots,1);
for iPlot = 1:nPlots
    plotObjs(iPlot) = plot(tArray,...
        sqrt(squeeze(intValsArray(iPlot,:))./max(intValsArray(:))),...
            '-','Color',colorList(iPlot,:),'LineWidth',1.5);
    hold on
    plot(tArray,sqrt(squeeze(intKinArray(iPlot,:))./max(intKinArray(:))),...
            '--','Color',colorList(iPlot,:),'LineWidth',1)
end
xlabel('Thickness (nm)')
ylabel('I^{1/2} (Arb)')
legend(plotObjs,peakNames)
xlim([0 40])

% Plot R factor between kinematical and multislice vs thickness

gam = 0.5;
Rarray = zeros(1,nUCs);
for iUC = 1:nUCs
    intsMS = intValsArray(:,iUC);
    intsKin = intKinArray(:,iUC);
    scaleFac = intsMS.^gam\intsKin.^gam;
    intsKin = (intsKin.^gam./scaleFac).^(1/gam);
    Rarray(iUC) = ...
        sum(abs(intsKin.^gam - intsMS.^gam))./sum(intsMS.^gam);
end

figure('Position',[100 100 300 200])
plot(tArray,100*Rarray,'k-',...
    'LineWidth',1.5,'MarkerSize',14)
xlabel('Thickness (nm)')
ylabel('R_{MS-Kin} (%)')
xlim([0 40])

save(['DiffAnalysis_' num2str(sDiff.E0/1e3) 'kV.mat'],...
    'tArray','Rarray','intValsArray','intKinArray','sDiff');

%% 3a. Compute real space electron waves
sDiff.storeRealSpaceCellEWs = true;
nUCs = 100; % Number of simulation cells to propagate through
theta_x = 0.0; % Incident angle of the e-beam (rad)
theta_y = 0.0;
coefs = [nUCs,0,theta_x,theta_y];
[~,EWImage,sDiff] = calcDiff(sDiff,coefs);

tArray = 0.1*sDiff.cellDim(3)*(1:nUCs); % Thickness (nm)
StackViewerEWs(abs(EWImage),tArray)

%% 3b. Visualize real space waves
% Show image of EW propagating through thickness
% ind1 = [384:-1:129 130:384];
% ind2 = [129*ones(1,256) 130:384];
% indx = sub2ind(size(EWImage),ind1,ind2);
% 
% EWlineProf = zeros(nUCs,numel(indx));
% for iUC = 1:nUCs
%     EWframe = EWImage(:,:,iUC);
%     EWlineProf(iUC,:) = EWframe(indx);
% end
% 
% figure;
% imagesc(sDiff.pot(1:512,1:512,1))

figure('Position',[100 100 350 400]);
imagesc([0 0.1*sDiff.cellDim(1)],[tArray(1) tArray(end)],...
    squeeze(abs(EWImage(64,:,:)))')
xlabel('x (nm)')
ylabel('z (nm)')
title('|\psi|')
colorbar
colormap(parula(1024))

figure('Position',[100 100 350 400]);
imagesc([0 0.1*sDiff.cellDim(1)],[tArray(1) tArray(end)],...
    squeeze(angle(EWImage(64,:,:)))')
xlabel('x (nm)')
ylabel('z (nm)')
title('\phi')
colorbar
colormap(hsv(1024))

%% 3c. Movie of real space envelope propagation

T = linspace(0,4,100);

[savefile,savepath] = uiputfile('*.avi');
saveName = [savepath savefile];

% create the video writer
generateMyColormaps
writerObj = VideoWriter(saveName);
writerObj.FrameRate = 10;

% open the video writer
open(writerObj);

scaling = max(squeeze(abs(EWImage(64,:,:))),[],[1 2]);
% write the frames to the video
figure('Position',[100 100 350 400]);
im = imagesc([0 0.1*sDiff.cellDim(1)],[tArray(1) tArray(end)],...
    squeeze(real(EWImage(64,:,:)))');
xlabel('x (nm)')
ylabel('z (nm)')
title('|\psi|')
colorbar
colormap(colormapBlWtRdPlus)
caxis([-scaling scaling])
for u=1:numel(T)
    
    set(im,'CData',...
        squeeze(real(EWImage(64,:,:).*exp(-2i*pi*T(u))))');
    writeVideo(writerObj, getframe(gcf));
end

% close the writer object
close(writerObj);
