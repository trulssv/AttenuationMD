function [a, aRef, seg, labels, pixelstomm] = LoadData(data_dir, phantom_case, phantom_type, visualization)

if nargin < 3
   visualization = 0; 
end


switch phantom_case
    case 1 % Small gammex I/Gd

        imgs_names = {'PE', 'PVC', 'I', 'Gd'};
        imgs = {'PE_img', 'PVC_img', 'm1_img', 'm2_img'};
        phantom_dir = 'gammex';
        data_paths = fullfile(data_dir, phantom_dir, 'phantomtype_1', imgs_names);
        BasisMaterials = {'polyethylene', 'pvc'};
        
        img_size = 512;
        recon_fov = 220;
        pixelstomm = recon_fov / img_size;

    case 2 % Small Gammex Ca/Gd
        imgs_filenames = {'PE', 'PVC', 'Ca', 'Gd'};
        imgs = {'PE_img', 'PVC_img', 'm1_img', 'm2_img'};
        phantom_dir = 'gammex';
        data_paths = fullfile(data_dir, phantom_dir, 'phantomtype_1', imgs_filenames);
        BasisMaterials = {'polyethylene', 'pvc'};
        
        img_size = 512;
        recon_fov = 220;
        pixelstomm = recon_fov / img_size;

    case 3 % Large gammex I/Gd
        imgs_filenames = {'PE', 'PVC', 'I', 'Gd'};
        imgs = {'PE_img', 'PVC_img', 'm1_img', 'm2_img'};
        phantom_dir = 'gammex';
        data_paths = fullfile(data_dir, phantom_dir, 'phantomtype_1', imgs_filenames);
        BasisMaterials = {'polyethylene', 'pvc'};
        
        img_size = 512;
        recon_fov = 220;
        pixelstomm = recon_fov / img_size;
        
    case 4 % Antropomorphic phantom 1
        imgs = {'PE', 'PVC', 'I', 'Gd'};
        imgs_ref = {'PE_ref', 'PVC_ref'};
        imgs_ref_names = {'PE_ref', 'PVC_ref'};
        
        phantom_dir = 'antropomorphic_phantoms';
        
        data_paths =fullfile(data_dir, phantom_dir, 'phantomtype_4', imgs);
        data_paths_ref = fullfile(data_dir, phantom_dir, 'phantomtype_4', imgs_ref);
        % Segmentation
        
        load(fullfile(data_dir, phantom_dir, 'phantomtype_4', 'segmentation'), 'seg');
        load(fullfile(data_dir, phantom_dir, 'phantomtype_4', 'labels'), 'names');
        labels = names;
        
        BasisMaterials = {'polyethylene', 'pvc'};
        cfg = GetCFG(phantom_case);
        img_size = 512;
        pixelstomm = cfg.recon_fov / img_size;
    case 5  
        imgs = {'PE', 'PVC', 'I', 'Gd'};
        imgs_ref = {'PE_ref', 'PVC_ref'};
        imgs_ref_names = {'PE_ref', 'PVC_ref'};

        phantom_dir = 'antropomorphic_phantoms';
        
        data_paths =fullfile(data_dir, phantom_dir, 'phantomtype_5', imgs);
        data_paths_ref = fullfile(data_dir, phantom_dir, 'phantomtype_5', imgs_ref);
        
        % segmentation
        
        load(fullfile(data_dir, phantom_dir, 'phantomtype_5', 'segmentation'), 'seg');
        load(fullfile(data_dir, phantom_dir, 'phantomtype_5', 'labels.mat'), 'names');
        
        labels = cellstr(names);
        
        BasisMaterials = {'polyethylene', 'pvc'};
        cfg = GetCFG(phantom_case);
        img_size = 512;
        pixelstomm = cfg.recon_fov / img_size;
    case 9
        imgs = {'PE', 'PVC', 'I', 'Gd'};
        imgs_ref = {'PE_ref', 'PVC_ref'};
        
        phantom_dir = 'antropomorphic_phantoms';
        
        data_paths =fullfile(data_dir, phantom_dir, 'phantomtype_9', imgs);
        data_paths_ref = fullfile(data_dir, phantom_dir, 'phantomtype_9', imgs_ref);
        
        % segmentation
        
        load(fullfile(data_dir, phantom_dir, 'phantomtype_4', 'segmentation'), 'seg');
        load(fullfile(data_dir, phantom_dir, 'phantomtype_4', 'labels.mat'), 'names');
        
        labels = cellstr(names);
        
        BasisMaterials = {'polyethylene', 'pvc'};
        cfg = GetCFG(phantom_case);
        img_size = 512;
        pixelstomm = cfg.recon_fov / img_size;
    case 50
        imgs = {'PE', 'PVC', 'I', 'Gd'};
        imgs_names = {'PE_norm', 'PVC_norm', 'Inorm', 'Gdnorm'};
        imgs_ref = {'PE_ref', 'PVC_ref'};
        imgs_ref_names = {'PE_ref', 'PVC_ref'};
        phantom_dir = 'antropomorphic_phantoms';
        phantom_case = "case_50";
        data_paths =fullfile(data_dir, phantom_dir, phantom_case, imgs);
        data_paths_ref = fullfile(data_dir, phantom_dir, phantom_case, imgs_ref);
        
        % segmentation
        
        load(fullfile(data_dir, phantom_dir, phantom_case, 'segmentation'), 'seg');
        load(fullfile(data_dir, phantom_dir, phantom_case, 'labels.mat'), 'names');
        
        labels = cellstr(names);
        
        BasisMaterials = {'polyethylene', 'pvc'};
        cfg = GetCFG(5);
        img_size = 512;
        pixelstomm = cfg.recon_fov / img_size;
end

%% Load Data

n_imgs = numel(imgs);

a = {};         % Basis images with contrast
aRef = {};      % Basis images without contrast

for i=1:n_imgs
   m = load(data_paths{i});
   m = m.(imgs_names{i});
   if i > 4
      m = 0.1 * m; 
   end
   a = [a, m];
end

switch phantom_type
    case 'gammex'
        aRef = [];
        seg = [];
        labels = [];
    case 'antropomorphic'

        n_imgs_ref = numel(imgs_ref);

        for i=1:n_imgs_ref
           m = load(data_paths_ref{i});
           m = m.(imgs_ref_names{i});
           aRef = [aRef, m];
        end

        %% Optional: Plot mono image

        if visualization
            plot_monos(a, aRef);
        end
end


