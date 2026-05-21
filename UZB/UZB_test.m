clear all, close all, clc

addpath ./../../adf_tools

%% Set path to data

dataroot = '/mnt/data/MDdata/UZB_megamonos/'; 
scan = {'CT99.449.2.2.1765289969.251640.', '303.12253.offline_input.'};


datapath_sino = strcat(dataroot, scan{1}, scan{2}, 'sino.adf');
datapath_img = strcat(dataroot, scan{1}, scan{2}, 'images.adf');
datapath_recon = strcat(dataroot, scan{1}, "DOCK_PIXELS.adf");


%recon = extract_data(dataroot, scan, 'recon');
%sino = extract_data(dataroot, scan, 'sino');
img = extract_data(dataroot, scan, 'mega_mono');

%% Extract slices

sz = size(sino);

nMats = sz(1);
nSlices = sz(3);

mid_slice = floor(nSlices / 2);

get_slice = @(sino, m, z) squeeze(sino(m, :, z, :));

material_sinos = struct;

for m=1:nMats
    m_sino = get_slice(sino, m, mid_slice);
    m_label = sprintf('material%d', m);
    material_sinos.(m_label) = m_sino;
end


%% plot data

plot_slices(material_sinos)

function plot_slices(data)
    
    labels = fieldnames(data);
    nMats = numel(labels);
    
    for m=1:nMats
        subplot(nMats, 1, m)
        
        label_m = labels{m};
        data_m = data.(label_m);

        imshow(data_m, [])
        h = gca; h.Visible = 'On';
        colorbar
        title(label_m)
        axis normal

    end
    

end





