%% Compute and display Kinematical tilt-averaged library
% as used in D.B. Durham et al, ArXiv Preprint 2022 

sDiff = setupSimKin();

nTheta = 64;
sigmaThetaMax = 0.16; % rad
sigmaThetaSamp = (sigmaThetaMax/nTheta)*(1:nTheta);

nUC = 100;
nIter = 10;

Ilib = computeTiltAveragedDiffraction(sigmaThetaSamp,nUC,nIter,sDiff);

tArray = 0.1*sDiff.cellDim(3)*(1:nUC);

StackViewerDiff(fftshift(fftshift(Ilib(:,:,:,end,5),1),2),tArray)

[savefile,savepath] = uiputfile('*.mat');
save([savepath savefile],'Ilib','sDiff',...
    'nTheta','sigmaThetaMax','sigmaThetaSamp',...
    'nUC','nIter','tArray');

%% Visualization

VisualizeTiltAveragedDiffraction;


