function MPctE = computeMPctEStack(IArray)
%COMPUTERSTACK Compute R for a stack of diffracted intensities referenced
%to the end of the stack
%   IArray - MxNxP array of peak intensities (order,thickness,parameter)

MPctE = 100*squeeze(...
    mean(abs(IArray(:,:,1:end) - IArray(:,:,end)),1) ...
    ./ sum(IArray(:,:,end),1));

end

