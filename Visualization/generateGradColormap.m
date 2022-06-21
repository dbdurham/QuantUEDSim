function cmap = generateGradColormap(cpts, N)
% GENERATEGRADCOLORMAP generates a gradient colormap using the array of
% color points defined in cpts
% cpts is an M x 4 array, containing R, G, B, and position (all in range 0
% to 1) for each of the M color points.
% N is the number of colors

% define a position vector representing the colors
p = linspace(0,1,N);
% linearly interpolate the colors between the points
cmap = zeros(N,3);
cmap(:,1) = interp1(cpts(:,4),cpts(:,1),p);
cmap(:,2) = interp1(cpts(:,4),cpts(:,2),p);
cmap(:,3) = interp1(cpts(:,4),cpts(:,3),p);
end
