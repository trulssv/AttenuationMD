% data_dir = './data/'; % Endoleak 559
data_dir = './data/mega_mono/CT99.449.2.2.1765289969.251640.303.12253'; % UZB I/Gd calibration insert
% data_dir = './data/mega_mono/CT99.456.2.1.1765364800.382548.302.12277'; % Mats 5 inserts with I/Gd
% data_dir = './data/mega_mono/CT99.560.2.1.1771925082.978325.303.13299'; %   Endoleak 560 artial
% data_dir = './data/mega_mono/CT99.561.2.1.1771925824.677157.303.13438'; % Endoleak 561
% data_dir = './data/mega_mono/CT99.563.2.1.1771927183.968347.303.13469'; % Endoleak 563
% data_dir = './data/mega_mono/CT99.564.2.1.1771927922.879175.304.13472'; % Endoleak 564
% data_dir = './data/mega_mono/CT99.565.2.1.1771928660.946249.303.13473';
% data_dir = './data/mega_mono/CT99.567.2.1.1771929984.766851.303.13295';


% ---- Load predicted maps ----
i_vol  = double(niftiread(fullfile(data_dir,'iodine_pred_vol.nii')));
gd_vol = double(niftiread(fullfile(data_dir,'gadolinium_pred_vol.nii')));

% Convert mass fractions to mg/mL
i_vol  = i_vol  * 4933;
gd_vol = gd_vol * 7900;

%% Define windowing parameters for display – cap values to avoid outliers dominating the color scaling

i_cap_max = 20;  % Cap predicted iodine at 20 mg/mL to avoid outliers dominating the color scaling
i_cap_min = 0;  % Allow some negative values for iodine to show in the color scaling (e.g., due to noise)
gd_cap_max = 10;  % Cap predicted gadolinium at 10 mg/mL to avoid outliers dominating the color scaling
gd_cap_min = 0;  % Allow some negative values for gadolinium to show in the color scaling (e.g., due to noise)
wa_cap_max = 1400;  % Cap predicted water at 2000 mg/mL (2000 kg/m3) to avoid outliers dominating the color scaling
wa_cap_min = 900;  % Allow some negative values for water to show in the color scaling (e.g., due to noise)
% Clamp negatives, cap iodine at 15 mg/mL
i_vol  = min(max(i_vol,  i_cap_min), i_cap_max);
gd_vol = min(max(gd_vol, gd_cap_min), gd_cap_max);

i_min  = i_cap_min;  i_max  = i_cap_max;
gd_min = gd_cap_min;  gd_max = gd_cap_max;

% ---- Load basis pair ----
ib_vol = double(niftiread(fullfile(data_dir,'iodine_basis_vol.nii')));
wb_vol = double(niftiread(fullfile(data_dir,'water_basis_vol.nii')));

% Convert basis maps to mg/mL
ib_vol = ib_vol * 4933;   % Iodine density 4933 kg/m3 -> mg/mL
wb_vol = wb_vol * 1000;   % Water density 1000 kg/m3  -> mg/mL

% Clamp negatives, cap basis iodine at 15 mg/mL
ib_vol = min(max(ib_vol, i_cap_min), i_cap_max);
wb_vol = min(max(wb_vol, wa_cap_min), wa_cap_max);

% Basis iodine shares the same display range as predicted iodine
ib_min = i_cap_min;  ib_max = i_cap_max;
wb_min = wa_cap_min;  wb_max = wa_cap_max;

% ---- Load GT maps if present ----
i_gt_file  = fullfile(data_dir, 'iodine_gt.nii');
gd_gt_file = fullfile(data_dir, 'gadolinium_gt.nii');
has_gt = isfile(i_gt_file) && isfile(gd_gt_file);

% ---- Load mask ----
mask_vol = double(niftiread(fullfile(data_dir, 'mask_vol.nii')));

% Build label volume: respect existing labels or derive via connected components
unique_labels = unique(mask_vol(:));
unique_labels(unique_labels == 0) = [];

if numel(unique_labels) > 1
    label_vol = mask_vol;
    roi_ids   = unique_labels;
else
    CC        = bwconncomp(mask_vol > 0, 6);
    label_vol = zeros(size(mask_vol));
    for k = 1:CC.NumObjects
        label_vol(CC.PixelIdxList{k}) = k;
    end
    roi_ids = (1:CC.NumObjects)';
end

% ---- Load GT maps ---------------------------------------------------
if ~has_gt
    error('No ground-truth files found in %s.\nRun label_ground_truth.m first.', data_dir);
end
i_gt_vol  = double(niftiread(i_gt_file));   % already in mg/mL
gd_gt_vol = double(niftiread(gd_gt_file));

% ---- Reload all volumes unclamped for statistics -------
% i_vol / gd_vol above were clamped to [0, cap]. Clamping negatives to 0
% biases the mean upward in low-signal ROIs (GT ≈ 0). Reload raw here so
% the statistical mean reflects the true prediction distribution.
i_raw  = double(niftiread(fullfile(data_dir, 'iodine_pred_vol.nii')))    * 4933;
gd_raw = double(niftiread(fullfile(data_dir, 'gadolinium_pred_vol.nii'))) * 7900;
wb_raw = double(niftiread(fullfile(data_dir, 'water_basis_vol.nii')))    * 1000;  % mg/mL, no clamping
ib_raw = double(niftiread(fullfile(data_dir, 'iodine_basis_vol.nii')))   * 4933;
% Voxels with water-basis below this threshold are classified as background.
wb_bg_threshold = 1000;   % mg/mL
ib_bg_threshold = 1;
%% ---- Middle axial slice --------------------------------------------
% Data layout is Z x X x Y  (dimension 1 = axial / Z)
nz  = size(i_vol, 1);
mid = round(nz / 2);

i_slice     = squeeze(i_vol(mid,    :, :));   % clamped – display only
gd_slice    = squeeze(gd_vol(mid,   :, :));   % clamped – display only
i_raw_slice  = squeeze(i_raw(mid,   :, :));   % unclamped – used for statistics
gd_raw_slice = squeeze(gd_raw(mid,  :, :));   % unclamped – used for statistics
wb_slice    = squeeze(wb_raw(mid,   :, :));
ib_slice    = squeeze(ib_raw(mid,   :, :));
i_gt_slice  = squeeze(i_gt_vol(mid, :, :));
gd_gt_slice = squeeze(gd_gt_vol(mid,:, :));
lbl_slice   = squeeze(label_vol(mid,:, :));

%% ---- Per-ROI statistics --------------------------------------------
n_rois       = numel(roi_ids);
i_gt_roi     = nan(n_rois, 1);
gd_gt_roi    = nan(n_rois, 1);
i_pred_roi   = nan(n_rois, 1);
gd_pred_roi  = nan(n_rois, 1);

for k = 1:n_rois
    roi_mask = (lbl_slice == roi_ids(k));

    if ~any(roi_mask(:))
        fprintf('ROI %d: no voxels on middle slice (z=%d) — skipped.\n', roi_ids(k), mid);
        continue
    end

    % GT is spatially constant within each ROI — extract from the first voxel
    i_gt_roi(k)  = i_gt_slice(find(roi_mask, 1));
    gd_gt_roi(k) = gd_gt_slice(find(roi_mask, 1));

    % Predicted: exclude background/air voxels using unclamped water-basis threshold
    tissue_mask = roi_mask & (wb_slice >= wb_bg_threshold | ib_slice >= ib_bg_threshold);
    if ~any(tissue_mask(:))
        fprintf('ROI %d: all voxels below water threshold (%.0f mg/mL) — using full ROI.\n', ...
                roi_ids(k), wb_bg_threshold);
        tissue_mask = roi_mask;
    end

    % Use unclamped values so that negative noise voxels are not forced to 0
    i_pred_roi(k)  = mean(i_raw_slice(tissue_mask));
    gd_pred_roi(k) = mean(gd_raw_slice(tissue_mask));

    % Diagnostics
    n_roi_vox    = nnz(roi_mask);
    n_tissue_vox = nnz(tissue_mask);
    fprintf('  ROI %d: %d vox in mask, %d pass threshold | I_raw: min=%.2f  med=%.2f  mean=%.2f  max=%.2f | Gd_raw: mean=%.2f\n', ...
            roi_ids(k), n_roi_vox, n_tissue_vox, ...
            min(i_raw_slice(tissue_mask)), median(i_raw_slice(tissue_mask)), ...
            i_pred_roi(k), max(i_raw_slice(tissue_mask)), gd_pred_roi(k));
end

% Keep only ROIs that had voxels on the middle slice
valid        = ~isnan(i_gt_roi);
i_gt_roi     = i_gt_roi(valid);
gd_gt_roi    = gd_gt_roi(valid);
i_pred_roi   = i_pred_roi(valid);
gd_pred_roi  = gd_pred_roi(valid);
roi_labels   = roi_ids(valid);
n_valid      = numel(roi_labels);

% Print summary table
fprintf('\n--- ROI summary (middle slice z=%d) ---\n', mid);
fprintf('%-6s  %-12s  %-12s  %-12s  %-12s\n', 'ROI', 'I_GT', 'I_pred', 'Gd_GT', 'Gd_pred');
for k = 1:n_valid
    fprintf('%-6d  %-12.3f  %-12.3f  %-12.3f  %-12.3f\n', ...
            roi_labels(k), i_gt_roi(k), i_pred_roi(k), gd_gt_roi(k), gd_pred_roi(k));
end

%% ---- Scatter plots -------------------------------------------------
i_lim  = [0,  max([i_gt_roi;  i_pred_roi;  eps]) * 1.15];
gd_lim = [0,  max([gd_gt_roi; gd_pred_roi; eps]) * 1.15];

figure('Name', 'Qualitative analysis: Predicted vs GT', ...
       'Units', 'normalized', 'Position', [0.10 0.15 0.70 0.60]);

% --- Iodine ---
subplot(1, 2, 1);
scatter(i_gt_roi, i_pred_roi, 80, 'b', 'filled'); hold on;
plot(i_lim, i_lim, 'k--', 'LineWidth', 1.2);   % identity line
for k = 1:n_valid
    text(i_gt_roi(k), i_pred_roi(k), sprintf('  ROI %d', roi_labels(k)), ...
         'FontSize', 8, 'Color', [0 0 0.7]);
end
xlim(i_lim); ylim(i_lim); axis square; grid on; box on;
xlabel('GT Iodine (mg/mL)');
ylabel('Predicted mean Iodine (mg/mL)');
title('Iodine');

% --- Gadolinium ---
subplot(1, 2, 2);
scatter(gd_gt_roi, gd_pred_roi, 80, 'r', 'filled'); hold on;
plot(gd_lim, gd_lim, 'k--', 'LineWidth', 1.2);
for k = 1:n_valid
    text(gd_gt_roi(k), gd_pred_roi(k), sprintf('  ROI %d', roi_labels(k)), ...
         'FontSize', 8, 'Color', [0.7 0 0]);
end
xlim(gd_lim); ylim(gd_lim); axis square; grid on; box on;
xlabel('GT Gadolinium (mg/mL)');
ylabel('Predicted mean Gadolinium (mg/mL)');
title('Gadolinium');

sgtitle(sprintf('Predicted vs GT — middle slice (z = %d) — %s', mid, data_dir), ...
        'Interpreter', 'none');

%% ---- ROI map on middle slice ---------------------------------------
% Apply display transform to match viewdata.m orientation
ref_display = rot90(fliplr(i_slice));
lbl_display = rot90(fliplr(lbl_slice));

ref_norm = mat2gray(ref_display, [i_cap_min, i_cap_max]);
rgb = repmat(ref_norm, [1 1 3]);

% Draw ROI borders in white
for k = 1:n_rois
    roi_mask_d = (lbl_display == roi_ids(k));
    if ~any(roi_mask_d(:)), continue; end
    border = bwperim(roi_mask_d);
    for ch = 1:3
        chan = rgb(:,:,ch);
        chan(border) = 1;
        rgb(:,:,ch) = chan;
    end
end

figure('Name', 'ROI map — middle slice', ...
       'Units', 'normalized', 'Position', [0.10 0.10 0.55 0.75]);
imshow(rgb, 'Border', 'tight'); hold on;

% Label each ROI at its centroid (computed in display space)
for k = 1:n_rois
    roi_mask_d = (lbl_display == roi_ids(k));
    if ~any(roi_mask_d(:)), continue; end
    props = regionprops(roi_mask_d, 'Centroid');
    cx = props(1).Centroid(1);
    cy = props(1).Centroid(2);
    text(cx, cy, num2str(roi_ids(k)), ...
         'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
title(sprintf('ROI map — middle slice (z = %d)', mid));
xlabel(data_dir, 'Interpreter', 'none');

%% ---- Per-ROI error plot --------------------------------------------
i_err  = i_pred_roi  - i_gt_roi;
gd_err = gd_pred_roi - gd_gt_roi;

x_ticks = 1:n_valid;
tick_labels = arrayfun(@(r) sprintf('ROI %d', r), roi_labels, 'UniformOutput', false);

figure('Name', 'Per-ROI prediction error', ...
       'Units', 'normalized', 'Position', [0.10 0.10 0.70 0.55]);

subplot(1, 2, 1);
bar(x_ticks, i_err, 'b'); hold on;
yline(0, 'k--', 'LineWidth', 1.2);
set(gca, 'XTick', x_ticks, 'XTickLabel', tick_labels);
ylabel('Predicted \minus GT (mg/mL)');
title('Iodine error'); grid on;

subplot(1, 2, 2);
bar(x_ticks, gd_err, 'r'); hold on;
yline(0, 'k--', 'LineWidth', 1.2);
set(gca, 'XTick', x_ticks, 'XTickLabel', tick_labels);
ylabel('Predicted \minus GT (mg/mL)');
title('Gadolinium error'); grid on;

sgtitle(sprintf('Per-ROI prediction error — middle slice (z = %d)', mid), 'Interpreter', 'none');
