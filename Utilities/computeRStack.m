function R = computeRStack(IArray)
%COMPUTERSTACK Compute R for a stack of diffracted intensities referenced
%to the end of the stack
%   IArray - MxNxP array of peak intensities (order,thickness,parameter)

Np = size(IArray,3);
Nt = size(IArray,2);
R = zeros(Nt,Np);
% alphaStore = zeros(Nt,Np);

for ip = 1:Np-1
    for it = 1:Nt
        I = IArray(:,it,ip);
        Iref = IArray(:,it,end);
        % Correct for scaling
        alpha = I\Iref;
%         alphaStore(it,ip) = alpha;
        I = I*alpha;
        % Compute R
        R(it,ip) = 100*sum(abs(I-Iref)) ./ sum(Iref);
    end
end


