function [sliceThickness] = computeSlices(zAtoms,zCell,tMin,numMax)
%COMPUTESLICES Compute slice thicknesses based on spacing between atoms 
%   zAtoms -- z positions of the atoms in the sim cell
%   zCell -- z height of the sim cell
%   tMin -- minimum thickness of a slice
%   numMax -- maximum number of slices to use

% handle inputs
if nargin < 4
    numMax = Inf;
end
if nargin < 3
    tMin = 0;
end

% 1. Identify the unique z positions of atoms
zAtoms = unique(zAtoms); % sorted in ascending order

% 2. Compute slice thicknesses
sliceThickness = [diff(zAtoms); zCell-sum(diff(zAtoms))];

% 3. Merge too-thin slices with neighboring slices

% while any slices are still too thin or there are too many left
while any(sliceThickness < tMin) | numel(sliceThickness) > numMax
    % find thinnest slice
    [~,iMin] = min(sliceThickness);
    % find thinnest neighboring slice
    if iMin == numel(sliceThickness)
        iMerge = iMin - 1;
    elseif iMin == 1
        iMerge = iMin + 1;
    else
        [~,iMerge] = max(sliceThickness([iMin-1,iMin+1]));
        iMerge = iMin + (iMerge*2)-3;
    end
    % add their values
    sliceThickness(iMerge) = sliceThickness(iMin)+sliceThickness(iMerge);
    % Remove thinnest slice
    sliceThickness = sliceThickness(~((1:numel(sliceThickness)) == iMin));
end



end

