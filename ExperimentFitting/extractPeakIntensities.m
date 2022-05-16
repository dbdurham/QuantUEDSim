function [intVals] = extractPeakIntensities(EW,qxa,qya,qxyPeaks,qUnit)
%EXTRACTPEAKINTENSITIES Extract peak intensities from simulated diff
%patterns
%   EW = M x N (complex) exit wave
%   qxa = M x N array of qx coordinates
%   qya = M x N array of qy coordinates
%   qxyPeaks = P x 2 array of peak positions (y, x)
%   qUnit = 'px','invAng'

if nargin < 2
    qUnit = 'px';
end

imageSize = size(EW);

nPeaks = size(qxyPeaks,1);
indPeaks = zeros(nPeaks,1); % linear indices
if strcmp(qUnit,'invAng') % 1/d (inverse Angstroms)
    for iPeak = 1:nPeaks
        [~,indPeaks(iPeak)] = ...
            min((qxyPeaks(iPeak,2)-qxa).^2 ...
            + (qxyPeaks(iPeak,1)-qya).^2);
    end
elseif strcmp(qUnit,'px') % 2D pixel indices
    indPeaks = sub2ind(imageSize,qxyPeaks(:,1),qxyPeaks(:,2));
end

% get simulated peak intensities (if defined)
intVals = abs(EW(indPeaks)).^2;

% % Damp values for thermal vibration
% damp = exp(-B*sRefine.exp_qxy(:,3));
% intVals(:) = intVals .* damp;

end

