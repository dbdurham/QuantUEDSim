%% Experimental h,k,l,I data for Au

[expfile,exppath] = uigetfile('*.mat');
load([exppath,expfile]);

setNameList = {'peakSet_200','peakSet_220','peakSet_400',...
    'peakSet_420','peakSet_440','peakSet_600','peakSet_620'};
peakNames = {'200','220','400','420','440','600','620'};
nSets = numel(setNameList);

hklExp = [2 0 0; 2 2 0; 4 0 0 ; 4 2 0; 4 4 0; 6 0 0; 6 2 0];
intsExp = zeros(nSets,1);
for iSet = 1:nSets
    intsExp(iSet) = ...
        mean(AnalysisDataOff.(setNameList{iSet}).volumeArrayNorm,[1 2]);
end

% Compute reciprocal space vectors
[~,~,~,uvwInit] = wyckoffGold();
Gvec = inv(uvwInit');
[GhklExp,GmagExp,~] = computeScatteringVectors(hklExp,Gvec);

% Show the avg off intensities
figure;
scatter(hklExp(:,1),hklExp(:,2),intsExp.^0.5,...
    intsExp.^0.5,'Filled','MarkerEdgeColor','k');
colormap(parula(1024))
colorbar()
set(gca,'ydir','reverse')
axis equal off
title('Measured static pattern I^{-1/2}')

save([exppath expfile(1:end-4) '_Ihkl.mat'],'peakNames','GmagExp','GhklExp','hklExp','intsExp')