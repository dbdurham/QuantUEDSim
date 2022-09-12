%% Compute R(%) vs thickness between kinematical and multislice 
%% across beam energies
% as used in D.B. Durham et al, ArXiv Preprint 2022 

%% Flat crystal

E0ListFlat = [30 40 50 60 70 80 95 ...
    110 130 150 175 200 300 400 575 750 ...
    1250 2000 3000 4000].*1e3;
nE = numel(E0ListFlat);

IArrayMSList = cell(nE,1);
IArrayKinList = cell(nE,1);
for iE = 1:nE
    options.E0 = E0ListFlat(iE);
    sDiffMS = setupSimMS(options);
    sDiffKin = setupSimKin(options);

    % Inputs for calculation
    theta1 = 0; % rad, x component of tilt
    theta2 = 0; % rad, y component of tilt
    nUC = 100; % Number of unit cells to simulate

    IDiffMS = calcDiffMS(theta1,theta2,nUC,sDiffMS);
    IDiffKin = calcDiffKin(-theta1,-theta2,nUC,sDiffKin);
    % Apply Debye-Waller Factor
    GmagStore = sqrt(sDiffKin.qxaStore.^2 + sDiffKin.qyaStore.^2);
    [~,DWFInt,~] = computeDWF(sDiffMS.uRMS,1,GmagStore);
    IDiffKin = IDiffKin.*DWFInt;
    
    % Extract beam intensities

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

    GhklTest = computeScatteringVectors(hklTest,sDiffKin.Gvec);

    IArrayMS = extractIntsFromDP(IDiffMS,...
        sDiffMS.qxaStore,sDiffMS.qyaStore,GhklTest);
    I0ArrayMS = extractIntsFromDP(IDiffMS,...
        sDiffMS.qxaStore,sDiffMS.qyaStore,[0 0 0]);
    IArrayKin = extractIntsFromDP(IDiffKin,...
        sDiffKin.qxaStore,sDiffKin.qyaStore,GhklTest);
    I0ArrayKin = extractIntsFromDP(IDiffKin,...
        sDiffKin.qxaStore,sDiffKin.qyaStore,[0 0 0]);

    tArray = (1:nUC)*sDiffMS.cellDim(3)*0.1; %nm

    showIvtPlusKin(IArrayMS,I0ArrayMS,...
        IArrayKin,tArray,peakNames);
    
    IArrayMSList{iE} = IArrayMS;
    IArrayKinList{iE} = IArrayKin;
end

%% Compute and display R vs thickness
figure('Position',[100 100 450 200]);

Rthresh = 10;
tThreshFlat = zeros(nE,1);
RindivFlat = zeros(nE,nUC,nPeaks);
tThreshIndivFlat = zeros(nE,nPeaks);
for iE = 1:nE
    
    IArrayMS = IArrayMSList{iE};
    IArrayKin = IArrayKinList{iE};
    Rflat = computeRStack(cat(3,IArrayKin,IArrayMS));
    semilogy(tArray,Rflat(:,1),'-','LineWidth',1.5)
    hold on
    xlabel('Thickness (nm)')
    ylabel('R (%)')
    
    tThreshFlat(iE) = tArray(find(Rflat>Rthresh,1));

    RindivFlat(iE,:,:) = 100*(abs(sqrt(IArrayKin)-sqrt(IArrayMS))...
        ./sqrt(IArrayMS))';
    for iPeak = 1:nPeaks
        tThreshIndivFlat(iE,iPeak) = ...
            tArray(find(RindivFlat(iE,:,iPeak)>Rthresh,1));
    end
    
end

hold on
plot([0 40],[10 10],'k--')

legendStr = {'30keV','40keV','50keV','60keV','70keV','80keV',...
    '95keV','110keV','130keV','150keV','175keV','200keV',...
    '300keV','400keV','575keV','750keV',...
    '1.25MeV','2MeV','3MeV','4MeV'};
legend(legendStr);
xlim([0 40])

colorList = jet(nPeaks).*0.8;
figure;
for iPeak = 1:nPeaks
    semilogx(E0ListFlat./1e3,tThreshIndivFlat(:,iPeak),'.',...
        'Color',colorList(iPeak,:),'MarkerSize',10)
    hold on
end
xlabel('Kinetic energy (eV)')
ylabel('Thickness where R > 10% (nm)')
legend(peakNames);

%% Load and compute R profiles for tilt-averaged simulations

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


folderPath = '../../DiffractionSimPackageOutputs/TestOutputs/';
Enames = {'30keV','40keV','50keV','60keV','70keV','80keV',...
    '95keV','110keV','130keV','150keV','175keV','200keV',...
    '300keV','400keV','575keV','750keV',...
    '1p25MeV','2MeV','3MeV','4MeV'};
E0List = [30 40 50 60 70 80 95 ...
    110 130 150 175 200 300 400 575 750 ...
    1250 2000 3000 4000].*1e3;
KinFileList = cellfun(@(c) [folderPath c '_u0p0894_Kin_10.mat'],...
    Enames,'UniformOutput',false);
MSFileList = cellfun(@(c) [folderPath c '_u0p0894_MS_01b.mat'],...
    Enames,'UniformOutput',false);
nFiles = numel(KinFileList);

RarrayList = cell(nFiles,1);
RindivList = cell(nFiles,1);
for iFile = 1:nFiles
    load(MSFileList{iFile})
    sDiff2 = sDiff;
    Ilib2 = Ilib;
    iEnd2 = size(Ilib2,5);
    
    load(KinFileList{iFile})
    sDiff1 = sDiff;
    Ilib1 = Ilib;
    iEnd = size(Ilib1,5);
    
    GhklTest = computeScatteringVectors(hklTest,sDiff1.Gvec);

    tArray = (1:nUC).*0.1*sDiff1.cellDim(3);

    % Apply Debye-Waller Factor
    GmagStore = sqrt(sDiff1.qxaStore.^2 + sDiff1.qyaStore.^2);
    [~,DWFInt,~] = computeDWF(sDiff2.uRMS,1,GmagStore);
    Ilib1 = Ilib1.*DWFInt;
    
    Rarray = zeros(nUC,nTheta);
    Rindiv = zeros(nUC,nTheta,nPeaks);
    for iTheta = 1:nTheta
        IArraySim2 = extractIntsFromDP(Ilib2(:,:,:,iTheta,iEnd2),...
            sDiff2.qxaStore,sDiff2.qyaStore,GhklTest);
        IArraySim1 = extractIntsFromDP(Ilib1(:,:,:,iTheta,iEnd),...
            sDiff1.qxaStore,sDiff1.qyaStore,GhklTest);
        Rstack = computeRStack(cat(3,...
            IArraySim1,IArraySim2));
        Rarray(:,iTheta) = Rstack(:,1);
        for iPeak = 1:nPeaks
            Rstack = 100*(abs(sqrt(IArraySim1(iPeak,:))-sqrt(IArraySim2(iPeak,:)))...
                ./sqrt(IArraySim2(iPeak,:)))';
            Rindiv(:,iTheta,iPeak) = Rstack(:,1);
        end    
    end
    RindivList{iFile} = Rindiv;
    RarrayList{iFile} = Rarray;

    cmap = violetFire.^0.5;

    figure;
    imagesc(sigmaThetaSamp([1 end])*1e3,...
        tArray([1 end]),...
        Rarray);
    hold on
    [X,Y] = meshgrid(sigmaThetaSamp.*1e3,tArray);
    contour(X,Y,Rarray,[3,3],'w--','LineWidth',1)
    contour(X,Y,Rarray,[10,10],'w--','LineWidth',1)
    xlabel('\sigma_{\theta} (mrad)')
    ylabel('Thickness (nm)')
    title('R_{BW - Kin} (%)')
    colormap(cmap)
    colorbar()
    caxis([0 100])
    set(gca,'ydir','normal')
end

%% Plot R vs thickness for rippled foil

iThetaToShow = 40;

figure('Position',[100 100 450 200]);

for iE = 1:nFiles
    
    R = RarrayList{iE};
    semilogy(tArray,R(:,iThetaToShow),'-','LineWidth',1.5)
    hold on
    xlabel('Thickness (nm)')
    ylabel('R (%)')

end

hold on
plot([0 40],[10 10],'k--')

legendStr = Enames;
legend(legendStr);
xlim([0 40])
ylim([1 100])

%% Plot R>10% vs sigmaTheta

Rthresh = 10;

figure('Position',[100 100 350 300]);

for iE = 1:nFiles
    R = RarrayList{iE};
    
    tThresh = zeros(1,nTheta);
    for iTheta = 1:nTheta
        iThresh = find(R(:,iTheta)>10,1);
        tThresh(iTheta) = tArray(iThresh);
    end
    
    plot(sigmaThetaSamp*1e3,tThresh,'-','LineWidth',1.5)
    hold on
    xlabel('\sigma_{\theta} (mrad)')
    ylabel('t_{R>10%} (nm)')

end

legendStr = Enames;
legend(legendStr);
xlim([0 160])
ylim([0 8])

%% Plot R>10% vs Energy for flat crystal, 100 mrad

Rthresh = 10;

iShowList = [10 20 40];
nShow = numel(iShowList);

figure('Position',[100 100 420 330]);

rippledColorList = [0.7 0 0;...
    0.7 0.7 0;...
    0.6 0.85 0.5];

tThresh = zeros(nFiles,nShow);

for iShow = 1:nShow
    iThetaToShow = iShowList(iShow);
    for iE = 1:nFiles
        R = RarrayList{iE};
        tThresh(iE,iShow) = tArray(find(R(:,iThetaToShow)>10,1));
    end
end

lobjs = gobjects(4,1);
lobjs(1) = semilogx(E0ListFlat./1e3,tThreshFlat,'ks','MarkerSize',6.5,...
    'MarkerFaceColor','k');

markerSymbolList = {'d','^','o'};
hold on
for iShow = 1:nShow
    lobjs(iShow+1) = semilogx(E0List./1e3,tThresh(:,iShow),...
        markerSymbolList{iShow},'MarkerSize',6,...
        'Color',rippledColorList(iShow,:),...
        'MarkerFaceColor',rippledColorList(iShow,:));
end
xlabel('Kinetic energy (keV)')
ylabel('Thickness where R_{Kin-Dyn} > 10% (nm)')
legendStr = arrayfun(@(c) ['\sigma_{\theta} = ' num2str(c) ' mrad'],...
    [0 sigmaThetaSamp(iShowList)*1e3],'UniformOutput',false);

legend(lobjs(end:-1:1),legendStr{end:-1:1},'NumColumns',1)

ax = gca;
set(ax,'Position',[0.1 0.2 0.84 0.7])

xlim([1e1 1e4])
ylim([0 8])

%% Plot R>10% for one sigma_theta, broken down by order

Rthresh = 10;

iThetaToShow = 40;

figure('Position',[100 100 420 300]);

colorList = jet(nPeaks).*0.8;

tThresh = zeros(nFiles,nPeaks);

for iPeak = 1:nPeaks
    for iE = 1:nFiles
        R = RindivList{iE};
        tThresh(iE,iPeak) = ...
            tArray(find(R(:,iThetaToShow,iPeak)>10,1));
    end
end

for iPeak = 1:nPeaks
    lobjs(iPeak) = semilogx(E0List./1e3,tThresh(:,iPeak),...
        '.','MarkerSize',10,...
        'Color',colorList(iPeak,:),...
        'MarkerFaceColor',colorList(iPeak,:));
    hold on
end
xlabel('Kinetic energy (keV)')
ylabel('Thickness where R_{Kin-Dyn} > 10% (nm)')

legend(lobjs,peakNames,'NumColumns',1)

ax = gca;
set(ax,'Position',[0.1 0.2 0.64 0.7])

xlim([5e1 1e4])
ylim([0 15])