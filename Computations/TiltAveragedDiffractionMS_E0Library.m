%% Compute and display Multislice tilt-averaged library
% as used in D.B. Durham et al, ArXiv Preprint 2022 

E0List = [30 40 50 60 70 80 95 110 130 175 200 300 400 575 ...
    750 1250 2000 3000 4000].*1e3;

[savefile,savepath] = uiputfile('*.mat');

options = struct;
options.uRMS = 0.0894; % 1D rms atomic displacement
options.qRange = 4; % q range to store in *output* (A^-1)
options.cellMult = 1;
options.downSampFacCell = 2;

nTheta = 64;
sigmaThetaMax = 0.16; % rad
sigmaThetaSamp = (sigmaThetaMax/nTheta)*(1:nTheta);

nUC = 100;
nIter = 9;

for E0 = E0List
    options.E0 = E0;
    sDiff = setupSimMS(options);

    % Compute the first set of tilt-averaged patterns
    Ilib = computeTiltAveragedDiffraction(sigmaThetaSamp,nUC,nIter,sDiff);
    % Further converge the smallest tilt-range patterns
    Ilib = computeTiltAveragedDiffraction(sigmaThetaSamp,nUC,10,sDiff,Ilib,10,4);

    tArray = 0.1*sDiff.cellDim(3)*(1:nUC);

    StackViewerDiff(fftshift(fftshift(Ilib(:,:,:,end,5),1),2),tArray)

    save([savepath num2str(E0./1e3) 'keV_' savefile],'Ilib','sDiff',...
        'nTheta','sigmaThetaMax','sigmaThetaSamp',...
        'nUC','nIter','tArray');
end

%% Visualization
VisualizeTiltAveragedDiffraction;


