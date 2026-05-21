clear all, close all, clc

%% Set up paths

addpath ./bowtie
addpath ./optim
addpath ./optim/utils
addpath ./segmentation
addpath ./utils
addpath ./visualization

%% Note: You need CATSIM for this to work!

catsim_path = '/media/data/catsim/Git_20250113_downloaded_2025_01_23/Git_20250113/main_GE';
addpath(genpath(catsim_path))

%% Set root data directory

data_dir = './data';
phantom_case = 50;   % only data for 1 and 50 in this directory 

%% Params

visualization = 1; % generate plots in ./visualization/plots
phantom_type = 'antropomorphic'; % gammex or antropomorphic

%% Load data

[a, aRef, seg, labels, pixelstomm] = LoadData(data_dir, phantom_case, phantom_type, visualization);
aC = {a{3}, a{4}};

%% Process data

[AB, AR, AC, MaskTot, detector_cfg] = PreprocessData(a, aRef, seg, labels, pixelstomm, phantom_type, visualization);

%% Solve 

AC_est = Solve(AB, AR);                   % Estimated contrast sinograms

if visualization
   plot_estimated_sinos(AC_est, AC)
end

%% Reconstruct

aC_est = Recon(AC_est, MaskTot, detector_cfg);          % Estimated contrast images     

if visualization
   plot_estimated_images(aC_est, aC)
end

