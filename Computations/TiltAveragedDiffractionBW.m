%% Compute and display Bloch Wave tilt-averaged library
% as used in D.B. Durham et al, ArXiv Preprint 2022 

options = struct;
options.E0 = 750e3; % Electron energy (eV)
options.uRMS = 0.0894; % 1D rms atomic displacement
options.GxyThresh = 4.5;
options.sThresh = 0.1;
options.cellMult = 2;
options.downSampFacCell = 2;
sDiff = setupSimBW(options);

nTheta = 64;
sigmaThetaMax = 0.16; % rad
sigmaThetaSamp = (sigmaThetaMax/nTheta)*(1:nTheta);

nUC = 100;
nIter = 9;

[savefile,savepath] = uiputfile('*.mat');

% Compute the first set of tilt-averaged patterns
Ilib = computeTiltAveragedDiffraction(sigmaThetaSamp,nUC,nIter,sDiff);
% Further converge the smallest tilt-range patterns
Ilib = computeTiltAveragedDiffraction(sigmaThetaSamp,nUC,nIter+1,sDiff,...
    Ilib,nIter+1,4);

tArray = 0.1*sDiff.cellDim(3)*(1:nUC);

StackViewerDiff(fftshift(fftshift(Ilib(:,:,:,end,5),1),2),tArray)

save([savepath savefile],'Ilib','sDiff',...
    'nTheta','sigmaThetaMax','sigmaThetaSamp',...
    'nUC','nIter','tArray');

%% Visualize tilt averaged diffraction information
VisualizeTiltAveragedDiffraction;