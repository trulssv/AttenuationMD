clear all, close all, clc

addpath(genpath('./../'))
addpath(genpath('./../../imtools3d'))

%% Set root data directory

data_dir = './data/mega_mono/';

dataroot = '/mnt/data/MDdata/UZB_megamonos/'; 

%% List all adf files

files = dir(fullfile(dataroot, '*.adf')); % List all files
N = numel(files);


for p=1:N
        fprintf('Processing file %d/%d...', p, N)
        datapath = strcat(dataroot, files(p).name);
        
        % Extract scan

        scan = split(files(p).name, '.');
        scan = join(scan(1:7), '.');
        scan = scan{1};
        scan_dir = strcat(data_dir, scan, '/');
        mkdir(scan_dir)
        
        data = double(extract_data(datapath, 'mega_mono'));

        io_unit = 0.02;
        wa_unit = 1;
        io_rho = 4933;
        wa_rho = 1000;
    
        data(2, :, :, :) = data(2, :, :, :) * io_unit / io_rho;
        data(1, :, :, :) = data(1, :, :, :) * wa_unit / wa_rho;
    

        N = size(data, 2);
        sz = size(data, 3);
        
        mask_vol = zeros([N, sz, sz]);
        i_basis_vol = zeros([N, sz, sz]);
        wa_basis_vol = zeros([N, sz, sz]);
        
        i_pred_vol = zeros([N, sz, sz]);
        gd_pred_vol = zeros([N, sz, sz]);
        mask = segment_slice_circular();
        for i=1:N
            fprintf('Processing slice %d / %d...', i, N)
            a_slice = {squeeze(data(1, i, :, :)), squeeze(data(2, i, :, :))};
            % plot_basis_maps(a, mask, i, N)
            aC = AttenuationMD(a_slice, mask);
            
            mask_vol(i, :, :) = mask;
            i_basis_vol(i, :, :) = a_slice{2};
            wa_basis_vol(i, :, :) = a_slice{1};
        
            i_pred_vol(i, :, :) = aC{1};
            gd_pred_vol(i, :, :) = aC{2};
        end
        
        
        niftiwrite(wa_basis_vol, strcat(scan_dir, 'water_basis_vol'))
        niftiwrite(i_basis_vol, strcat(scan_dir, 'iodine_basis_vol'))
        niftiwrite(mask_vol, strcat(scan_dir, 'mask_vol'))
        niftiwrite(i_pred_vol, strcat(scan_dir, 'iodine_pred_vol'))
        niftiwrite(gd_pred_vol, strcat(scan_dir, 'gadolinium_pred_vol'))
end



function mask = segment_slice_circular()
    
    %% Filter out background and air
    centers = [
        422.5000 157.8333;
        290.5000 61.8333;
        297.1864 133.3608;
        344.3537 169.0063;
        225.8952 149.2032;
        354.7954 245.6983;
        235.2567 223.0148;
        279.5436 255.7799;
        290.7053 335.7124;
        160.3650 235.6167;
        127.9599 279.9037;
        170.4466 314.4691;
        216.1737 343.9937;
        259.7405 379.9993;
        93.7546 322.7504;
        137.6814 354.4353;
        185.6915 390.8546;
        228.1766 421.7199;
        ];

    P = size(centers, 1);
    sample_radius = 22;
    radii = ones([P, 1]) * sample_radius;

    mask = circles2mask(centers, radii, [512, 512]);
    

    

    
end

function mask = circles2mask(centers, radii, imSize)
% circles2mask Create binary mask from multiple circles
%
% Inputs:
%   centers : Nx2 array [x, y]
%   radii   : Nx1 array
%   imSize  : [rows, cols]
%
% Output:
%   mask    : logical image

    % Image grid
    [X, Y] = meshgrid(1:imSize(2), 1:imSize(1));

    % Initialize mask
    mask = false(imSize);

    % Loop over circles
    for k = 1:size(centers,1)
        cx = centers(k,1);
        cy = centers(k,2);
        r  = radii(k);

        % Circle equation
        circle = (X - cx).^2 + (Y - cy).^2 <= r^2;

        % Combine
        mask = mask | circle;
    end
end
