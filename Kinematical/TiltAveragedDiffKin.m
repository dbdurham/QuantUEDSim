%% Calculate tilt-averaged kinematical pattern over range of thicknesses

%% Calculate initial kinematical pattern

calcDiffKinematical_01_DBD;

%% Calculate tilt-averaged patterns

UCmax = 100;
dCrystal = (1:1:UCmax).*cellDim(3);
nUC = length(dCrystal);
% Beam characteristics
E0 = 3000; %keV
lamElec = wavelengthElectrons(E0); % Electron de Broglie wavelength (Angstroms)
% Other inputs
Vcryst = prod(cellDim); % Unit cell volume in Angstroms^3
alpha = 0; % Inverse absorption depth (Angstroms^-1)
B = 0; % Debye-Waller B Factor

% NOTE: in principle more efficient to calculate extinction distances in
% advance... but computation is fast enough here

% Could symmetrize DPs, but presently set up for calculating arbitrary hkl

% % define tilt range
sigmaThetaMax = 160e-3;
sigmaThetaSamp = (2.5:2.5:160) * 1e-3;

func = @(x,y) computeWeightedIntsKin(...
    x,y,... % sample tilt angles
    sigmaThetaSamp,sigmaThetaMax,... % tilt spreads to study
    dCrystal,... % thicknesses
    Ghkl,SFMag,... % g vectors and structure factors of reflections
    lamElec,... % electron wavelength
    Vcryst,... % unit cell volume
    alpha,... % absorption length
    B,... % Debye-Waller Factor
    false); % Symmetrize DPs in x and y

% Initialize libraries
nTheta = numel(sigmaThetaSamp);
nIter = 9;
Ilib = zeros(nOrders,...
    nUC,...
    nTheta,...
    nIter);
epsLib = zeros(nUC,nTheta,nIter);

iStart = 1;
iEnd = nIter;

for iIter = iStart:iEnd
    tic
    if iIter == 1
        Ilib(:,:,:,1) = extendTrapz2D(func,-3*sigmaThetaMax,3*sigmaThetaMax,...
            -3*sigmaThetaMax,3*sigmaThetaMax,Ilib(:,:,:,1),iIter);
    else
        Ilib(:,:,:,iIter) = extendTrapz2D(func,-3*sigmaThetaMax,3*sigmaThetaMax,...
            -3*sigmaThetaMax,3*sigmaThetaMax,Ilib(:,:,:,iIter-1),iIter);
    end
    
    if iIter > 1
        % Fractional change as function of thickness and tilt spread
        epsLib(:,:,iIter) = squeeze(sum(abs(...
            Ilib(:,:,:,iIter) - Ilib(:,:,:,iIter-1)),1)...
            ./ sum(Ilib(:,:,:,iIter-1),1));
    end
        
    disp(['Computed iteration #' num2str(iIter)])
    toc
    disp(['Fractional change in diff stack at max spread: ' num2str(mean(epsLib(:,end,iIter)))])
end

[savefile,savepath] = uiputfile('*.mat');
save([savepath savefile],'sigmaThetaMax','sigmaThetaSamp','UCmax',...
    'Ilib','epsLib','iEnd','nIter',...
    'dCrystal','hkl','Ghkl','SFMag','lamElec','Vcryst','alpha','B','cellDim')

%% Compute convergence

nTheta = numel(sigmaThetaSamp);
nUC = length(dCrystal);


epsLib = zeros(nUC,nTheta,nIter);
dThetaList = 0.75*sigmaThetaMax./2.^(0:nIter-1);

useSimp = false;

for iIter = 3:nIter
    IstorePrev = Ilib(:,:,:,iIter-1);
    Istore = Ilib(:,:,:,iIter);
    
    if useSimp
        IstorePrev2 = Ilib(:,:,:,iIter-2);
        % Apply Simpson's Rule
        Istore = (16/15)*Istore - (1/15)*IstorePrev;
        IstorePrev = (16/15)*IstorePrev - (1/15)*IstorePrev2;
    end
    
    clear('Istack2');
    % Fractional change as function of thickness and tilt spread
    epsLib(:,:,iIter) = squeeze(sum(abs(Istore - IstorePrev),1)...
        ./ sum(IstorePrev,1));
end

tArray = 0.1*cellDim(3)*(1:UCmax);
iToShow = 3:nIter;
nToShow = numel(iToShow);
cmap = jet(nToShow).*0.8;
iTheta = 8;

figure;
for iIter = iToShow
    semilogy(tArray,squeeze(epsLib(:,iTheta,iIter)),'LineWidth',1.5,'Color',cmap(iIter+1-iToShow(1),:))
    hold on
end
xlabel('Thickness (nm)')
ylabel('Fractional change')
title(['Thickness dependent change for \sigma_{\theta} = ' ...
    num2str(sigmaThetaSamp(iTheta)*1e3) 'mrad'])
legend(arrayfun(@(c) ['d\theta = ' num2str(c*1e3) 'mrad'],dThetaList(iToShow),...
    'UniformOutput',false))
xlim([0 105])

%% 2D plot of intensities

nUC = 42;
tUC = 0.1*cellDim(3)*nUC;
iTheta = 32;
I = Ilib(:,nUC,iTheta,end);

figure;
scatter(Ghkl(:,1),Ghkl(:,2),50,I,'filled','MarkerEdgeColor','k')
xlabel('q_{x} (Angstroms^{-1})')
ylabel('q_{y} (Angstroms^{-1})')
title(['Tilt-averaged pattern for \sigma_{\theta} = ' ...
    num2str(sigmaThetaSamp(iTheta)*1e3) ' mrad, t = ' num2str(tUC,3) ' nm'])
colormap('hot')
colorbar
Isort = sort(I(:));
caxis([0 Isort(end-1)])

%% Select diffracted beams to study

% List of q values of peaks to study
h1 = [0.4902 0]; 
h2 = [0 0.4902];
qxyPeaks = [1 0; ...
    1 1; 2 0; ...
    2 1; 2 2; ...
    3 0; 3 1; ...
    3 2]...
    *[h1; h2];
peakNames = {'200',...
    '220','400',...
    '420','440',...
    '600','620',...
    '640'};

nPeaks = size(qxyPeaks,1);
indPeaks = zeros(nPeaks,1);
for iPeak = 1:nPeaks
    [~,indPeaks(iPeak)] = ...
        min((qxyPeaks(iPeak,1)-Ghkl(:,1)).^2 ...
        + (qxyPeaks(iPeak,2)-Ghkl(:,2)).^2);
end

%% Plot intensity with increasing thickness

nUC = 50;

I0Array = zeros(1,nUC);
intValsThickness = zeros(nPeaks,nUC);
tArray = (1:nUC)*cellDim(3)*0.1; %nm

iTheta = 1;

Istack = Ilib(:,:,iTheta,end);

for iUC = 1:nUC
    I = Istack(:,iUC);
    intValsThickness(:,iUC) = I(indPeaks);
end
intValsNorm = intValsThickness ./ max(intValsThickness,[],2);
intValsRel = intValsThickness ./ max(intValsThickness,[],1);

% Plot the peak intensities as function of thickness
normMode = 'rel'; 
peaksToShow = 1:nPeaks;

switch normMode
    case 'raw'
        y = intValsThickness(peaksToShow,:);
        labelY = 'Intensity (Arb)';
    case 'norm'
        y = intValsNorm(peaksToShow,:);
        labelY = 'Normalized Intensity';
    case 'rel'
        y = intValsRel(peaksToShow,:);
        labelY = 'Relative Intensity';
end

figure('Position',[200 200 450 200]);

plot(tArray,y,'.-')
xlabel('Thickness (nm)')
ylabel(labelY)
title(['Warped sample \sigma_{\theta} = ' ...
    num2str(sigmaThetaSamp(iTheta)*1e3) ' mrad'])
legend(peakNames{peaksToShow})

% figure('Position',[200 200 450 200]);
% plot(tArray,I0Array,'k.-')
% xlabel('Thickness (nm)')
% ylabel('Peak intensity')
% title('Warped sample \sigma_{\theta} = 30 mrad')