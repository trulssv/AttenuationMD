% label_ground_truth.m
%
% Interactive script to assign ground-truth I / Gd concentrations (mg/mL)
% to each cylindrical ROI in every processed scan under data/mega_mono/.
%
% Workflow
% --------
%   For each scan directory that contains a mask_vol.nii:
%     1. Load mask and a reference image (predicted Gd map).
%     2. For each labelled ROI (connected component) show its location
%        overlaid on the reference image.
%     3. Prompt the user to enter [I_conc, Gd_conc] for that ROI.
%     4. After all ROIs are labelled, write iodine_gt.nii and
%        gadolinium_gt.nii into the same scan directory.
%
%   If iodine_gt.nii and gadolinium_gt.nii already exist, they are loaded
%   so the user can review each ROI, edit the values, or press Enter to
%   keep the current values unchanged.

clear; close all; clc;

data_root = fullfile(fileparts(mfilename('fullpath')), 'data', 'mega_mono');

scan_dirs = dir(data_root);
scan_dirs = scan_dirs([scan_dirs.isdir] & ~startsWith({scan_dirs.name}, '.'));

for s = 1:numel(scan_dirs)

    scan_path = fullfile(data_root, scan_dirs(s).name);
    mask_file = fullfile(scan_path, 'mask_vol.nii');

    if ~isfile(mask_file)
        fprintf('Skipping %s (no mask_vol.nii)\n', scan_dirs(s).name);
        continue
    end

    gt_i_file  = fullfile(scan_path, 'iodine_gt.nii');
    gt_gd_file = fullfile(scan_path, 'gadolinium_gt.nii');

    fprintf('\n========================================\n');
    fprintf('Scan: %s\n', scan_dirs(s).name);
    fprintf('========================================\n');

    answer = input('Label this scan? [Y/n]: ', 's');
    if strcmpi(strtrim(answer), 'n')
        fprintf('Skipped.\n');
        continue
    end

    %% Load volumes
    mask_vol = double(niftiread(mask_file));   % Z x X x Y  (int or float)

    % Load a reference image for visual guidance (prefer iodine basis, fall back
    % to water basis)
    ref_candidates = {'iodine_basis_vol.nii', 'water_basis_vol.nii'};
    ref_vol = [];
    for c = 1:numel(ref_candidates)
        f = fullfile(scan_path, ref_candidates{c});
        if isfile(f)
            ref_vol = double(niftiread(f));
            ref_name = ref_candidates{c};
            break
        end
    end
    if isempty(ref_vol)
        error('No reference image found in %s', scan_path);
    end

    % Convert Gd / I pred from mass-fraction (stored units) to mg/mL for display
    if contains(ref_name, 'gadolinium')
        ref_vol = ref_vol * 7900;   % Gd density kg/m3 -> mg/mL
        ref_clim = [0 10];
    elseif contains(ref_name, 'iodine')
        ref_vol = ref_vol * 4933;   % I  density kg/m3 -> mg/mL
        ref_clim = [0 20];
    else
        ref_vol = ref_vol * 1000;   % Water
        ref_clim = [0 1200];
    end

    %% Find ROIs as connected components in the mask
    % mask_vol may be binary or labelled. Handle both cases.
    unique_labels = unique(mask_vol(:));
    unique_labels(unique_labels == 0) = [];   % remove background

    if isempty(unique_labels)
        fprintf('  Mask contains no non-zero voxels. Skipping.\n');
        continue
    end





    % Build a label volume (same size as mask) where each connected component
    % in any non-zero mask region gets a unique integer label.
    % If already labelled (>1 unique values), respect existing labels.
    if numel(unique_labels) > 1
        label_vol = mask_vol;               % already labelled
        roi_ids   = unique_labels;
    else
        % Binary mask: connected-component analysis to separate inserts
        CC = bwconncomp(mask_vol > 0, 6);   % 6-connectivity for 3-D
        label_vol = zeros(size(mask_vol));
        for k = 1:CC.NumObjects
            label_vol(CC.PixelIdxList{k}) = k;
        end
        roi_ids = (1:CC.NumObjects)';
    end

    n_rois = numel(roi_ids);
    fprintf('  Found %d ROI(s).\n', n_rois);

    %% Per-ROI labelling loop
    has_existing_gt = isfile(gt_i_file) && isfile(gt_gd_file);

    if has_existing_gt
        fprintf('  Existing ground-truth files found. Loading current ROI values.\n');
        iodine_gt     = double(niftiread(gt_i_file));
        gadolinium_gt = double(niftiread(gt_gd_file));
    else
        iodine_gt     = zeros(size(mask_vol));
        gadolinium_gt = zeros(size(mask_vol));
    end

    % Store per-ROI answers so we can go back and correct them.
    i_concs  = nan(n_rois, 1);
    gd_concs = nan(n_rois, 1);
    for k = 1:n_rois
        roi_mask_k = (label_vol == roi_ids(k));
        i_vals = iodine_gt(roi_mask_k);
        gd_vals = gadolinium_gt(roi_mask_k);

        if ~isempty(i_vals)
            i_concs(k) = median(i_vals(:), 'omitnan');
        end
        if ~isempty(gd_vals)
            gd_concs(k) = median(gd_vals(:), 'omitnan');
        end
    end

    fig = figure('Name', sprintf('Ground-truth labelling — %s', scan_dirs(s).name), ...
                 'Units', 'normalized', 'Position', [0.05 0.05 0.85 0.85]);

    k = 1;
    while k <= n_rois

        roi_mask = (label_vol == roi_ids(k));

        % Find the axial slice with the most ROI voxels for display
        slice_counts = squeeze(sum(sum(roi_mask, 2), 3));  % Z x 1
        [~, best_slice] = max(slice_counts);

        ref_slice  = squeeze(ref_vol(best_slice, :, :));
        roi_slice  = squeeze(roi_mask(best_slice, :, :));

        % Apply the same display transform used in viewdata.m so the
        % orientation matches the ground-truth overlay there.
        ref_slice = rot90(fliplr(ref_slice));
        roi_slice = rot90(fliplr(roi_slice));

        %% Visualise
        clf(fig);
        subplot(1, 2, 1);
        imagesc(ref_slice, ref_clim);
        colormap(gca, 'hot'); colorbar; axis image off;
        title(sprintf('%s – slice %d', strrep(ref_name,'_',' '), best_slice));

        subplot(1, 2, 2);
        % Overlay: reference image + ROI border in cyan
        ref_norm = (ref_slice - ref_clim(1)) / (ref_clim(2) - ref_clim(1));
        ref_norm = repmat(mat2gray(ref_norm), [1 1 3]);
        overlay  = ref_norm;
        border   = bwperim(roi_slice);
        overlay(:,:,1) = overlay(:,:,1) - 0.4 * border;  % reduce red on border
        overlay(:,:,2) = overlay(:,:,2) + 0.4 * border;  % increase green
        overlay(:,:,3) = overlay(:,:,3) + 0.4 * border;  % increase blue → cyan
        overlay  = max(0, min(1, overlay));
        imshow(overlay); axis image off;
        title(sprintf('ROI %d of %d  (highlighted in cyan)', k, n_rois));
        sgtitle(scan_dirs(s).name, 'Interpreter', 'none');
        drawnow;

        % Show previously entered value as hint (if revisiting)
        if ~isnan(i_concs(k))
            hint = sprintf(' [current: %.2f  %.2f]', i_concs(k), gd_concs(k));
        else
            hint = '';
        end

        %% Ask for concentrations
        while true
            if k > 1
                back_hint = '  Type ''back'' to go back to the previous ROI.\n';
            else
                back_hint = '';
            end
            fprintf(back_hint);
            if ~isnan(i_concs(k)) && ~isnan(gd_concs(k))
                prompt = sprintf('  ROI %d/%d%s — enter [Iodine_mg_per_mL, Gd_mg_per_mL] or press Enter to keep unchanged: ', k, n_rois, hint);
            else
                prompt = sprintf('  ROI %d/%d%s — enter [Iodine_mg_per_mL, Gd_mg_per_mL] (e.g. 5 0): ', k, n_rois, hint);
            end
            raw = strtrim(input(prompt, 's'));

            if strcmpi(raw, 'back')
                if k > 1
                    k = k - 1;
                    fprintf('  Went back to ROI %d.\n', k);
                else
                    fprintf('  Already at the first ROI.\n');
                end
                break   % re-enter the while-k loop for the (new) k
            end

            if isempty(raw)
                if ~isnan(i_concs(k)) && ~isnan(gd_concs(k))
                    roi_mask_k = (label_vol == roi_ids(k));
                    iodine_gt(roi_mask_k)     = i_concs(k);
                    gadolinium_gt(roi_mask_k) = gd_concs(k);
                    fprintf('  ROI %d unchanged: I = %.2f mg/mL,  Gd = %.2f mg/mL\n', k, i_concs(k), gd_concs(k));
                    k = k + 1;
                    break
                end

                fprintf('  No existing value for this ROI. Please enter two non-negative numbers.\n');
                continue
            end

            vals = str2double(strsplit(raw));
            if numel(vals) == 2 && all(~isnan(vals)) && all(vals >= 0)
                i_concs(k)  = vals(1);
                gd_concs(k) = vals(2);

                roi_mask_k = (label_vol == roi_ids(k));
                iodine_gt(roi_mask_k)     = i_concs(k);
                gadolinium_gt(roi_mask_k) = gd_concs(k);

                fprintf('  ROI %d: I = %.2f mg/mL,  Gd = %.2f mg/mL\n', k, i_concs(k), gd_concs(k));
                k = k + 1;
                break
            end
            fprintf('  Invalid input. Please enter two non-negative numbers.\n');
        end

    end

    close(fig);

    %% Save ground-truth volumes
    % NIfTI info from the mask for consistent header
    mask_info = niftiinfo(mask_file);
    
    gt_info = mask_info;
    gt_info.Datatype = 'single';
    gt_info.BitsPerPixel = 32;

    niftiwrite(single(iodine_gt),  gt_i_file,  gt_info, 'Compressed', false);
    niftiwrite(single(gadolinium_gt), gt_gd_file, gt_info, 'Compressed', false);

    fprintf('  Saved iodine_gt.nii and gadolinium_gt.nii to\n  %s\n', scan_path);
end

fprintf('\nDone.\n');
