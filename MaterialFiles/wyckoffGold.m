function [atoms,cellDim,abcCell,uvwInit] = wyckoffGold(options)

% This script generates orthogonal output cells of FCC Au

% Inputs
flagPlot = true;

% Lattice vectors
aVec = 4.08; % angstroms (Quoted for Ted Pella single crystal foils)

% Atomic numbers 
Z = 79;

% x y z atomic_# class_# class_letter
% classes:
% [1 1]
sites = [0 0 0   Z 1 1;
    0 0.5 0.5    Z 1 1;
    0.5 0 0.5    Z 1 1;
    0.5 0.5 0    Z 1 1;];
% 7 free wyckoff positions in the space group

cellMult = 1; % Number of unit cells along one dim of sim cell

if nargin > 0
    if isfield(options,'sites') && ~isempty(options.sites)
        sites = options.sites;
    end

    if isfield(options,'cellMult') && ~isempty(options.cellMult)
        cellMult = options.cellMult;
    end
end

xProj = cellMult.*[1 0 0];
yProj = cellMult.*[0 1 0];
zProj = [0 0 1];
UCtile = [cellMult cellMult 1]; % Number of cells to tile along each dimension
% units of cellDim to fit all
abcShift = [1/4 1/4 0.001]; % Shift so all atoms are contained within the unit cell

% Build cell
numAtoms = sum(sites(:,5));
ID = zeros(numAtoms,1);
abcCell = zeros(numAtoms,3);
for a0 = 1:size(sites,1) % For each atomic site...
    inds = (1:sites(a0,5)) + sum(sites(1:(a0-1),5));
    
    ID(inds) = sites(a0,4);
    abcCell(inds,:) = sites(a0,1:3);
end
u0 = [aVec 0 0];
v0 = [0 aVec 0];
w0 = [0 0 aVec];
uvwInit = [u0;v0;w0];
% Projection vectors
u0 = uvwInit(1,:);
v0 = uvwInit(2,:);
w0 = uvwInit(3,:);
u = xProj(1)*u0 + xProj(2)*v0 + xProj(3)*w0;
v = yProj(1)*u0 + yProj(2)*v0 + yProj(3)*w0;
w = zProj(1)*u0 + zProj(2)*v0 + zProj(3)*w0;

% cell dimensions
cellDim = [norm(u) norm(v) norm(w)];
%
% tile sites
[ya,xa,za] = meshgrid( ...
    0:(UCtile(2)-1),0:(UCtile(1)-1),0:(UCtile(3)-1));
p = [xa(:) ya(:) za(:)];
[bInd,pInd] = meshgrid(1:size(abcCell,1),1:size(p,1));
IDtile = ID(bInd(:));
abcInit = abcCell(bInd(:),:) + p(pInd(:),:);
% Project sites into new unit cell
% initial cartesian coordinates
xyzInit = (uvwInit' * abcInit')';
abcProj = ([u;v;w]' \ xyzInit')';
abcProj(:) = mod(abcProj,1);
atoms = [ ...
    abcProj(:,1)*cellDim(1) ...
    abcProj(:,2)*cellDim(2) ...
    abcProj(:,3)*cellDim(3) ...
    IDtile];
%
for iDim = 1:3
    atoms(:,iDim) = ...
        mod(atoms(:,iDim)+abcShift(iDim)*cellDim(iDim),cellDim(iDim));
end

%
if flagPlot == true
    
    % plot output
    figure(1)
    clf
    set(gcf,'color','w')
    hold on
    line([0 0 cellDim(2) cellDim(2) 0],...
        [0 cellDim(1) cellDim(1) 0 0],...
        [0 0 0 0 0],...
        'linewidth',1,'color','k');
    line([0 0 cellDim(2) cellDim(2) 0],...
        [0 cellDim(1) cellDim(1) 0 0],...
        [0 0 0 0 0]+cellDim(3),...
        'linewidth',1,'color','k');
    
    s = atoms(:,4) == Z;
    scatter3(atoms(s,2),atoms(s,1),atoms(s,3),...
        'marker','o','sizedata',30,'linewidth',1,...
        'markerfacecolor','w','markeredgecolor','r')
    
%     figObjHand = findobj(gcf,'Children','axes_of_data');
% direction = [1 0 0];
% rotate(s,direction,25)
    hold off
    axis equal off
%     view([0 0 1])
%     view([0 0 0.])
    xlim([0 cellDim(2)]+[-1 1])
    ylim([0 cellDim(1)]+[-1 1])
    zlim([0 cellDim(3)]+[-1 1])
    
    
    % Cell angle
    gamma = acos(sum(u.*w) / norm(u) / norm(w));
    
    % cell volume
    volInit = abs(sum(cross(u0,v0).*w0));
    volNew = abs(sum(cross(u,v).*w));
    volRatio = volNew / volInit;
    
    disp(['gamma = ' sprintf('%.04f',gamma*180/pi) ' degrees'])
    disp(['volume ratio = ' sprintf('%.02f',volRatio)])
    
end
end


