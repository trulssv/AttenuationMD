function aC = Recon(AC, MaskTot, detector_cfg)



%% Reconstruct images from sinos


aC = {};

for i=1:numel(AC)
   sino = AC{i}; 
   img = 1/detector_cfg.pixelstomm * ifanbeam(...,
    sino, detector_cfg.D_px, ...
    'OutputSize', detector_cfg.outputsize, ...
    'FanSensorGeometry', 'arc', ...
    'FanSensorSpacing', detector_cfg.spacing, ...
    'FanRotationIncrement', detector_cfg.rot_increment, ...
    'Interpolation', 'linear', ...
    'Filter', 'Ram-Lak');
    aC = [aC, img];
end

Ipred_img = aC{1};
Gdpred_img = aC{2};


%% Post-process imgs

Ipred_img0 = Ipred_img .* MaskTot;
Gdpred_img0 = Gdpred_img .* MaskTot;

Ipred_img = zeros(size(Ipred_img0));
Gdpred_img = zeros(size(Gdpred_img0));

CC = bwconncomp(MaskTot);
numRegions = CC.NumObjects;

for k = 1:numRegions
    submask = false(size(MaskTot));
    submask(CC.PixelIdxList{k}) = true;

    Gdmask = Gdpred_img0 .* submask;
    Imask = Ipred_img0 .* submask;
    Ipred_img = Ipred_img + maskedLowPass(Imask, submask, 3);
    Gdpred_img = Gdpred_img + maskedLowPass(Gdmask, submask, 3);
end

aC = {Ipred_img, Gdpred_img};