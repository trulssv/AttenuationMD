clear all, close all, clc

%% Set up paths

addpath(genpath('./../'))

%% Set root data directory

data_dir = './data/';


%% Params

visualization = 0; % generate plots in ./visualization/plots

%% Load data

rho_w = 1;          % g / cm3
rho_i = 4.933;      % g / cm3

% water

wa = load(strcat(data_dir, 'water.mat'));  % presumably in mg/ml
wa = wa.i1 / 1000;    % g / ml = g / cm3
wa = wa / rho_w;

% iodine

io = load(strcat(data_dir, 'iodine.mat')); % presumably in 0.1mg/ml
io = io.i2 / 1e5;     % g / ml = g / cm3
io = io / rho_i;

pixelstomm = 200 / 512; % From dicominfo


%% Segmentation step

load(strcat(data_dir, 'mask.mat'), 'mask')


mask(wa<0.9) = 0; % air seg mask

wB = wa;
wB(mask) = 1;
wR = wa - wB;

iB = io;
iB(mask) = 0;
iR = io - iB;

aB = {wB, iB};
aR = {wR, iR};
a = [aB, aR];

%% Calculate sinograms

% GT sinogtrams

%% Fanbeam


detector_cfg = struct();

D = 625;
D_px = D / pixelstomm;
spacing = 0.1;
rot_increment = 1;
outputsize = size(a{1}, 1);

detector_cfg.pixelstomm = pixelstomm;
detector_cfg.D = 625;
detector_cfg.D_px = D / pixelstomm;
detector_cfg.spacing = 0.1;
detector_cfg.rot_increment = 1;
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

aC_est = Recon(AC_est, mask, detector_cfg);          % Estimated contrast images     

iEst = aC_est{1};
gdEst = aC_est{2};

subplot(2, 2, 1)
imshow(iEst, [])
colorbar
title('Iodine map (predicted)')


subplot(2, 2, 2)
imshow(gdEst, [])
colorbar
title('gadolinium map (predicted)')

subplot(2, 2, 3)
imshow(wa, [])
colorbar
title('Water map (initial)')


subplot(2, 2, 4)
imshow(io, [])
colorbar
title('Iodine map (initial)')


