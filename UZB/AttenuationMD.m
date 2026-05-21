function aC = AttenuationMD(a, mask)


%% Assuming the basis a is given by Water/Iodine



%% Calculate sinograms
wB = a{1};
wB(mask) = 1;
wR = a{1} - wB;

iB = a{2};
iB(mask) = 0;
iR = a{2} - iB;

aB = {wB, iB};
aR = {wR, iR};
a = [aB, aR];

% GT sinograms

%% Fanbeam geometry


detector_cfg = struct();

D = 625;
spacing = 0.1;
rot_increment = 1;
pixelstomm = 200 / 512; % From dicominfo

detector_cfg.pixelstomm = pixelstomm;
detector_cfg.D = D;
detector_cfg.D_px = D / pixelstomm;
detector_cfg.spacing = spacing;
detector_cfg.rot_increment = rot_increment;
detector_cfg.outputsize = size(a{1}, 1);

A = {};

for i=1:numel(a)
    [s, gamma, rotAngles] = fanbeam(a{i}, detector_cfg.D_px, ...,
    'FanSensorGeometry', 'arc', ...
    'FanSensorSpacing', detector_cfg.spacing, ...
    'FanRotationIncrement', detector_cfg.rot_increment);
    
    s = detector_cfg.pixelstomm * s;
    A = [A, s];
end


AB = A(1:2);
AR = A(3:4);

%% Solve 

AC_est = Solve(AB, AR);                   % Estimated contrast sinograms


%% Reconstruct

aC = Recon(AC_est, mask, detector_cfg);          % Estimated contrast images     