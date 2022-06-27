function sDiff = computePotential(sDiff)

correctForTilt = false;
if correctForTilt
    theta_x_corr = cos(sDiff.theta_x);
    theta_y_corr = cos(sDiff.theta_y);
else
    theta_x_corr = 1;
    theta_y_corr = 1;
end

% 1. Calculate projected potentials for each atom type 
% Prepare variables
xyLeng = ceil(sDiff.potBound./sDiff.pixelSize);
xvec = -xyLeng(1):xyLeng(1);
yvec = -xyLeng(2):xyLeng(2);
xr = xvec*sDiff.pixelSize(1)*theta_x_corr; 
yr = yvec*sDiff.pixelSize(2)*theta_y_corr;
% Lookup table for atom types
sDiff.atomTypes = unique(sDiff.atoms(:,4));
sDiff.potLookup = zeros( ...
    length(xvec),length(yvec),length(sDiff.atomTypes));
for a0 = 1:length(sDiff.atomTypes)
    sDiff.potLookup(:,:,a0) = projPot(sDiff.atomTypes(a0),xr,yr);
end
sDiff.potLookupFFT = fft2(sDiff.potLookup);
sDiff.xvec = xvec;
sDiff.yvec = yvec;
% Subpixel potential shift array
[qxa,qya] = makeFourierCoords(size(sDiff.potLookup),1);
sDiff.potShiftX = -2i*pi*qxa;
sDiff.potShiftY = -2i*pi*qya;
% Potential coords
qxPot = makeFourierCoords(size(sDiff.potLookup,1),...
    sDiff.pixelSize(1)*theta_x_corr);
qyPot = makeFourierCoords(size(sDiff.potLookup,2),...
    sDiff.pixelSize(2)*theta_y_corr);
[qyPot,qxPot] = meshgrid(qyPot,qxPot);
sDiff.q2Pot = qxPot.^2 + qyPot.^2;

% 2. Build the potential slices
sDiff.pot = ...
    zeros([sDiff.imageSize sDiff.numSlices sDiff.numFP]);
for a0 = 1:sDiff.numSlices
    inds = find(sDiff.sliceIndex == a0);

    for a1 = 1:sDiff.numFP % For each frozen lattice configuration

        for a2 = 1:length(inds) % For each atom
            [~,indType] = min(abs(sDiff.atomTypes ...
                - sDiff.atoms(inds(a2),4)));
            x = sDiff.atoms(inds(a2),1) / sDiff.pixelSize(1);
            y = sDiff.atoms(inds(a2),2) / sDiff.pixelSize(2);
            % positions in pixels

            % apply subpixel shifts
            xR = round(x);
            yR = round(y);
            if sDiff.numFP > 1 % Apply random shifts if using frozen phonon model
                dx = x - xR + randn*sDiff.uFP/sDiff.pixelSize(1);
                dy = y - yR + randn*sDiff.uFP/sDiff.pixelSize(2);
            else
                dx = x - xR;
                dy = y - yR;
            end

            % Write into output potential
            xp = mod(sDiff.xvec+xR,sDiff.imageSize(1))+1;
            yp = mod(sDiff.yvec+yR,sDiff.imageSize(2))+1;
            [Xp,Yp] = meshgrid(xp,yp);
            if sDiff.DWDamping
                [~,~,DWFAmp] = computeDWF(sDiff.uRMS,1,sqrt(sDiff.q2Pot));
                DWFilter = DWFAmp;
            else
                DWFilter = ones(size(sDiff.potLookupFFT,[1 2]));
            end
            potToAdd = ifft2(sDiff.potLookupFFT(:,:,indType) ... % FFT of atomic potential
                .* exp(dx*sDiff.potShiftX + dy*sDiff.potShiftY)... % subpixel shift
                .* DWFilter ... % Debye-Waller Damping
                ,'symmetric');

            sDiff.pot(:,:,a0,a1) = sDiff.pot(:,:,a0,a1) + ...
                accumarray([Xp(:) Yp(:)],potToAdd(:),sDiff.imageSize);

        end
    end

end

% 3. Generate and store the transmission function
% NOTE: when accounting for tilt, sigma should be multiplied by 
% k/kz = 1/cos(theta), but I have already accounted for this factor
% in the potentials by leaving out a cos(theta) correction there.
if sDiff.flagUseAntiAliasing && strcmp(sDiff.AntiAliasingMode,'2/3') % Band limit if needed
    sDiff.expPot = ifft2(fft2(exp(1i*sDiff.sigma*sDiff.pot)).*sDiff.qMask);
else
    sDiff.expPot = exp(1i*sDiff.sigma*sDiff.pot);
end

end