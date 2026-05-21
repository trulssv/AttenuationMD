clear all, close all, clc

addpath(genpath('./../'))
addpath(genpath('./../../imtools3d'))

%% Set root data directory

data_dir = './data/mega_mono/';

%% pipeline


% 1. specify directory
% 2. Load every dcm in relevant range
% 3. Generate masks (or manually if no other options exist)
% 4. Save I/Wa basis maps and ROI masks
% 5. Run MD algorithm on all files

data_format = 'mega_mono';




switch data_format
    case 'dicom'
        
        dataroot = '/mnt/data/MDdata/Images_Endoleak/20260224-e559/s310/';
        processdir = '/mnt/data/MDdata/Images_Endoleak/20260224-e559/s310_processed';

        files = dir(fullfile(dataroot, '*.dcm')); % List all files
        N = numel(files);
        sz = 512;
        
        mask_vol = zeros([N, sz, sz]);
        i_basis_vol = zeros([N, sz, sz]);
        wa_basis_vol = zeros([N, sz, sz]);
        
        i_pred_vol = zeros([N, sz, sz]);
        gd_pred_vol = zeros([N, sz, sz]);
        
        
        
        
        
        
        for i=1:N
            fprintf('Processing slice %d / %d...', i, N)
            datapath = strcat(dataroot, files(i).name);
            info = dicominfo(datapath);
            a = load_basis_maps(info);
            mask = segment_slice(a);
            % plot_basis_maps(a, mask, i, N)
            aC = AttenuationMD(a, mask);
            
            mask_vol(i, :, :) = mask;
            i_basis_vol(i, :, :) = a{2};
            wa_basis_vol(i, :, :) = a{1};
        
            i_pred_vol(i, :, :) = aC{1};
            gd_pred_vol(i, :, :) = aC{2};
        end
        
        
        niftiwrite(wa_basis_vol, strcat(data_dir, 'water_basis_vol'))
        niftiwrite(i_basis_vol, strcat(data_dir, 'iodine_basis_vol'))
        niftiwrite(mask_vol, strcat(data_dir, 'mask_vol'))
        niftiwrite(i_pred_vol, strcat(data_dir, 'iodine_pred_vol'))
        niftiwrite(gd_pred_vol, strcat(data_dir, 'gadolinium_pred_vol'))
    
    case 'mega_mono'
        dataroot = '/mnt/data/MDdata/UZB_megamonos/'; 
        % UZB I/Gd Cal
        scan = {'CT99.448.2.1.1765289335.300826.', '302.offline_input.'};

        %scan = {'CT99.449.2.2.1765289969.251640.', '303.12253.offline_input.'};

        %% Scan 
        % smaller phantom with 5 I/Gd samples
        % scan = {'CT99.456.2.1.1765364800.382548.', '302.12277.offline_input.'};
        

        % scan = {'CT99.560.2.1.1771925082.978325.', '303.13299.offline_input.'};  % Endoleak arterial
        % scan = {'CT99.561.2.1.1771925824.677157.', '303.13438.offline_input.'};  % 561 - Endoleak venous
        % scan = {'CT99.563.2.1.1771927183.968347.', '303.13469.offline_input.'};  % 563 - Endoleak I only
        % scan = {'CT99.564.2.1.1771927922.879175.', '304.13472.offline_input.'};   % 564 - Endoleak I only
        %scan = {'CT99.565.2.1.1771928660.946249.', '303.13473.offline_input.'};    % 565 - 
        %scan = {'CT99.567.2.1.1771929984.766851.', '303.13295.offline_input.'};
        
        s = split(scan{1}, '.');
        scan_id = str2num(s{2});


        scan_dir = join(scan, '');
        scan_dir = strcat(data_dir, scan_dir{1}, '/');

        datapath_img = strcat(dataroot, scan{1}, scan{2}, 'images.adf');
        
        data = double(extract_data(datapath_img, 'mega_mono'));

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
        mask = segment_slice_circular(scan_id);
        for i=1:N
            fprintf('Processing slice %d / %d...', i, N)
            a_slice = {squeeze(data(1, i, :, :)), squeeze(data(2, i, :, :))};
            mask = segment_slice(a_slice);
            % plot_basis_maps(a, mask, i, N)
            aC = AttenuationMD(a_slice, mask);
            
            mask_vol(i, :, :) = mask;
            i_basis_vol(i, :, :) = a_slice{2};
            wa_basis_vol(i, :, :) = a_slice{1};
        
            i_pred_vol(i, :, :) = aC{1};
            gd_pred_vol(i, :, :) = aC{2};
        end

        %% Save results
        
        mkdir(scan_dir)
        niftiwrite(wa_basis_vol, strcat(scan_dir, 'water_basis_vol'))
        niftiwrite(i_basis_vol, strcat(scan_dir, 'iodine_basis_vol'))
        niftiwrite(mask_vol, strcat(scan_dir, 'mask_vol'))
        niftiwrite(i_pred_vol, strcat(scan_dir, 'iodine_pred_vol'))
        niftiwrite(gd_pred_vol, strcat(scan_dir, 'gadolinium_pred_vol'))



end


function plot_basis_maps(a, mask, i, N)
    clf()
    subplot(2, 2, 1)
    imshow(a{1}, [0, 1.5])
    colorbar
    title(sprintf('water map: slice %d / %d', i, N))

    subplot(2, 2, 2)
    imshow(a{2}, [0, 0.002])
    colorbar
    title(sprintf('iodine map: slice %d / %d', i, N))

    subplot(2, 2, 3)
    imshow(mask)
    colorbar
    title(sprintf('ROI mask: slice %d / %d', i, N))
    shg

 
end


function a = load_basis_maps(dcm_info)

    wa_field = 'Private_0053_1095';
    io_field = 'Private_0053_1096';
    
    % preprocess

    wa = reshape(double(dcm_info.(wa_field)), [512, 512]);
    io = reshape(double(dcm_info.(io_field)), [512, 512]);

    % There is some bug in the dtype conversion, so certain negative
    % numbers end up as very large positive numbers. Therefore, we set all
    % values larger than som threshold to zero. 

    threshold = 5e3;
    
    wa(wa>threshold) = 0; % in mg/ml
    io(io>threshold) = 0; % in 0.02 mg/ml

    io_unit = 0.02;
    wa_unit = 1;
    io_rho = 4933;
    wa_rho = 1000;
    
    io = io * io_unit / io_rho;
    wa = wa * wa_unit / wa_rho;

    a = {wa, io};


    
end


function mask = segment_slice(a)
    
    %% Filter out background and air

    range_w = [0.94, 1.02];
    range_i = [2 * 5.5e-4, 2 * 8.5e-4];
    w_min = 0.94;
    i_min = 2 * 1e-4;

    m_w = (a{1} > range_w(1) & a{1} < range_w(2)) | a{1} < w_min;
    m_i = (a{2} > range_i(1) & a{2} < range_i(2)) | a{2} < i_min;

    m = ~(m_w & m_i);
    

   %% Only consider values inside phantom

   x_range = [80, 380];
   y_range = [50, 450];
    
  
    phantom_mask = false(size(m));
    
    phantom_mask(y_range(1):y_range(2), x_range(1):x_range(2)) = true;

    m = m & phantom_mask;

    %% Only keep large connected components

    cc = bwconncomp(m, 8);  % your binary image
    threshold = 250;      % minimum area
    
    mask = false(size(m));  % initialize output mask
    
    for i = 1:cc.NumObjects
        if numel(cc.PixelIdxList{i}) > threshold
            mask(cc.PixelIdxList{i}) = true;
        end
    end
    

    
end

function mask = segment_slice_circular(scan_id)
    if nargin == 0
        scan_id = 449;
    end
    %% Filter out background and air

    switch scan_id

        case 448


            centers = [
               195.5000 65.5000
               115.5000 139.0000
               190.5000 218.0000
               265 297.5000
               343.5000 222.5000
               299.5000 415.0000
            
            ];

            sample_radius = 20;

        case 449
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
            sample_radius = 22;
        case 456
            centers = [
               238.9842 188.4579;
               399.3211 219.4474;
               367.3211 380.4579;
               208.3316 351.1526;
               119.7421 225.8474;
            ];
            sample_radius = 22;
        case 560
            centers = [
               238.9842 188.4579;
               399.3211 219.4474;
               367.3211 380.4579;
               208.3316 351.1526;
               119.7421 225.8474;
            ];
            sample_radius = 22;

    end

    P = size(centers, 1);
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


% 456 
% e560 contrast arterial
% 561 venous arterial




