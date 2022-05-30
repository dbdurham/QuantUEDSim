function RBands = computeRBands(IArray,tArray,tBands)
%COMPUTERBANDS Compute R(%) for varying params
%   IArray = MxNxP array of peak intensities (order,thickness,parameter)
%   tArray = 1xN array of thickness values
%   tBands = Qx2 array of lower and upper bounds for desired thickness
%   bands

nTests = size(IArray,3);
nBands = size(tBands,1);

Rend = computeRStack(IArray);

RBands = zeros(nBands,nTests);
for iBand = 1:nBands
    indLow = find(tArray > tBands(iBand,1),1,'first');
    indHigh = find(tArray < tBands(iBand,2),1,'last');
    RBands(iBand,:) = mean(Rend(indLow:indHigh,:),1);
end

end

