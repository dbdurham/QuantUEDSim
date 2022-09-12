%% Compute temperature-dependent Bloch Wave tilt-averaged libraries
% as used in D.B. Durham et al, ArXiv Preprint 2022

options = struct;
options.E0 = 750e3; % Electron energy (eV)
options.GxyThresh = 4.5;
options.sThresh = 0.1;
options.cellMult = 2;
options.downSampFacCell = 2;

du2List = 0:0.001:0.015;
u2Init = 0.024;
uRMSList = sqrt((u2Init+du2List)/3);
nSims = numel(uRMSList);

nTheta = 64;
sigmaThetaMax = 0.16; % rad
sigmaThetaSamp = (sigmaThetaMax/nTheta)*(1:nTheta);

nUC = 50;
nIter = 8;
tArray = 0.1*sDiff.cellDim(3)*(1:nUC);

[savefile,savepath] = uiputfile('*.mat');

for iSim = 2:nSims
    options.uRMS = uRMSList(iSim);
    sDiff = setupSimBW(options);

    % Compute the first set of tilt-averaged patterns
    Ilib = computeTiltAveragedDiffraction(sigmaThetaSamp,nUC,nIter,sDiff);
    % Further converge the smallest tilt-range patterns
    Ilib = computeTiltAveragedDiffraction(sigmaThetaSamp,nUC,nIter+1,sDiff,...
        Ilib,nIter+1,4);

    save([savepath savefile(1:end-4) ...
        '_u' pad(num2str(iSim),2,'left','0') '.mat'],...
        'Ilib','sDiff',...
        'nTheta','sigmaThetaMax','sigmaThetaSamp',...
        'nUC','nIter','tArray');
end

%% Visualization
VisualizeTiltAveragedDiffraction;

