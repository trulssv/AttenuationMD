function viewdata(plane, mode)
%VIEWDATA  Visualise one plane of material-decomposition volumes.
%
%   MODE 'combined' (default) – RGB overlay, two panels:
%     Left panel  – Predicted maps:  Iodine (blue) + Gadolinium (red)
%     Right panel – Basis pair:      Iodine (blue) + Water (green)
%
%   MODE 'separate' – individual colormapped panels on a 2×2 grid:
%     Top-left:     Predicted Iodine      (blue colormap)
%     Bottom-left:  Predicted Gadolinium  (red colormap)
%     Top-right:    Basis Iodine          (blue colormap)
%     Bottom-right: Basis Water           (green colormap)
%
%   If iodine_gt.nii and gadolinium_gt.nii exist in the scan directory,
%   a compare button cycles the right-hand view between the basis pair,
%   ground-truth I/Gd, and absolute error maps |prediction - GT|. This
%   works in both combined and separate modes.
%
%   Both modes share the same colormaps, display ranges, and slice slider.
%   A "Switch to …" button in the top-right corner toggles between modes
%   at any time without reloading data or losing the current slice.
%
%   viewdata()                        – axial, combined
%   viewdata('axial')                 – axial, combined
%   viewdata('coronal', 'separate')   – coronal, separate
%   viewdata('sagittal','combined')   – sagittal, combined









if nargin < 1, plane = 'axial';    end
if nargin < 2, mode  = 'combined'; end

if ~ismember(lower(mode), {'combined','separate'})
    error('mode must be ''combined'' or ''separate''');
end

% data_dir = './data/'; % Endoleak 559
% data_dir = './data/mega_mono/CT99.449.2.2.1765289969.251640.303.12253'; % UZB I/Gd calibration insert
data_dir = './data/mega_mono/CT99.456.2.1.1765364800.382548.302.12277'; % Mats 5 inserts with I/Gd
% data_dir = './data/mega_mono/CT99.560.2.1.1771925082.978325.303.13299'; %   Endoleak 560 artial
% data_dir = './data/mega_mono/CT99.561.2.1.1771925824.677157.303.13438'; % Endoleak 561
% data_dir = './data/mega_mono/CT99.563.2.1.1771927183.968347.303.13469'; % Endoleak 563
% data_dir = './data/mega_mono/CT99.564.2.1.1771927922.879175.304.13472'; % Endoleak 564
% data_dir = './data/mega_mono/CT99.565.2.1.1771928660.946249.303.13473';
% data_dir = './data/mega_mono/CT99.567.2.1.1771929984.766851.303.13295';
data_dir = './data/mega_mono/CT99.448.2.1.1765289335.300826.302';

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
compare_view = 'basis';

if has_gt
    i_gt_vol  = double(niftiread(i_gt_file));
    gd_gt_vol = double(niftiread(gd_gt_file));
    i_gt_vol  = min(max(i_gt_vol,  i_cap_min), i_cap_max);
    gd_gt_vol = min(max(gd_gt_vol, gd_cap_min), gd_cap_max);
    i_err_vol  = min(abs(i_vol  - i_gt_vol),  i_cap_max);
    gd_err_vol = min(abs(gd_vol - gd_gt_vol), gd_cap_max);
else
    i_gt_vol  = [];
    gd_gt_vol = [];
    i_err_vol  = [];
    gd_err_vol = [];
end

% Data layout is Z x X x Y
[nz, nx, ny] = size(i_vol);

switch lower(plane)
    case 'axial',    n_slices = nz;  init = round(nz/2);  % dim 1 = Z
    case 'coronal',  n_slices = nx;  init = round(nx/2);  % dim 2 = X
    case 'sagittal', n_slices = ny;  init = round(ny/2);  % dim 3 = Y
    otherwise, error('plane must be ''axial'', ''coronal'', or ''sagittal''');
end

% ---- Colormaps ----
ncolors = 256;
ramp    = linspace(0, 1, ncolors)';
cmap_blue  = [zeros(ncolors,2), ramp];
cmap_red   = [ramp, zeros(ncolors,2)];
cmap_green = [zeros(ncolors,1), ramp, zeros(ncolors,1)];

% ---- Figure (created once; layout rebuilt on toggle) ----
fig = figure('Name', sprintf('%s Viewer', plane), ...
    'Units','normalized','Position',[0.05 0.1 0.9 0.82], ...
    'Color','k');

% Shared mutable state (written by buildLayout / toggleMode)
current_slice = init;
s   = [];
lbl = [];
% Axes handles – only the relevant set is valid for the current mode
ax_pred = [];  ax_basis = [];
ax_ip   = [];  ax_gp    = [];  ax_ib = [];  ax_wb = [];

buildLayout();
updateSlice();

% ======== nested functions ========

    % ------------------------------------------------------------------
    % buildLayout  –  (re)creates all axes, colorbars, slider and button.
    %                 Called on startup and every time mode is toggled.
    % ------------------------------------------------------------------
    function buildLayout()
        clf(fig);
        set(fig, 'Color', 'k', 'Name', sprintf('%s Viewer – %s', plane, getModeLabel()));

        if strcmpi(mode, 'combined')
            % ---- Two-panel RGB overlay layout ----
            ax_pred  = axes('Position',[0.03 0.10 0.39 0.82], 'Color','k', ...
                            'XColor','w','YColor','w');
            ax_basis = axes('Position',[0.52 0.10 0.39 0.82], 'Color','k', ...
                            'XColor','w','YColor','w');

            makeCB([0.43 0.52 0.022 0.37], cmap_blue,  i_min,  i_max,  'Iodine pred (mg/mL)');
            makeCB([0.43 0.10 0.022 0.37], cmap_red,   gd_min, gd_max, 'Gadolinium pred (mg/mL)');
            if strcmp(compare_view, 'gt')
                makeCB([0.91 0.52 0.022 0.37], cmap_blue,  i_min,  i_max,  'Iodine GT (mg/mL)');
                makeCB([0.91 0.10 0.022 0.37], cmap_red,   gd_min, gd_max, 'Gadolinium GT (mg/mL)');
            elseif strcmp(compare_view, 'error')
                makeCB([0.91 0.52 0.022 0.37], cmap_blue,  0, i_max,  'Iodine abs. error (mg/mL)');
                makeCB([0.91 0.10 0.022 0.37], cmap_red,   0, gd_max, 'Gadolinium abs. error (mg/mL)');
            else
                makeCB([0.91 0.52 0.022 0.37], cmap_blue,  ib_min, ib_max, 'Iodine basis (mg/mL)');
                makeCB([0.91 0.10 0.022 0.37], cmap_green, wb_min, wb_max, 'Water basis (mg/mL)');
            end

        else  % separate – left: predicted (I/Gd), right: basis (I/Wa)
            ax_ip = axes('Position',[0.03 0.53 0.37 0.38], 'Color','k', ...
                         'XColor','w','YColor','w');
            ax_gp = axes('Position',[0.03 0.10 0.37 0.38], 'Color','k', ...
                         'XColor','w','YColor','w');
            ax_ib = axes('Position',[0.51 0.53 0.37 0.38], 'Color','k', ...
                         'XColor','w','YColor','w');
            ax_wb = axes('Position',[0.51 0.10 0.37 0.38], 'Color','k', ...
                         'XColor','w','YColor','w');

            makeCB([0.41 0.53 0.022 0.38], cmap_blue,  i_min,  i_max,  'Iodine pred (mg/mL)');
            makeCB([0.41 0.10 0.022 0.38], cmap_red,   gd_min, gd_max, 'Gadolinium pred (mg/mL)');
            if strcmp(compare_view, 'gt')
                makeCB([0.89 0.53 0.022 0.38], cmap_blue,  i_min,  i_max,  'Iodine GT (mg/mL)');
                makeCB([0.89 0.10 0.022 0.38], cmap_red,   gd_min, gd_max, 'Gadolinium GT (mg/mL)');
            elseif strcmp(compare_view, 'error')
                makeCB([0.89 0.53 0.022 0.38], cmap_blue,  0, i_max,  'Iodine abs. error (mg/mL)');
                makeCB([0.89 0.10 0.022 0.38], cmap_red,   0, gd_max, 'Gadolinium abs. error (mg/mL)');
            else
                makeCB([0.89 0.53 0.022 0.38], cmap_blue,  ib_min, ib_max, 'Iodine basis (mg/mL)');
                makeCB([0.89 0.10 0.022 0.38], cmap_green, wb_min, wb_max, 'Water basis (mg/mL)');
            end
        end

        % ---- Toggle buttons (top-right) ----
        uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.68 0.955 0.14 0.036], 'String', 'Switch mode', ...
            'BackgroundColor',[0.25 0.25 0.25], 'ForegroundColor','w', ...
            'FontSize',9, 'Callback', @(~,~) toggleMode());

        if has_gt
            if strcmp(compare_view, 'basis')
                compare_str = 'Show GT';
            elseif strcmp(compare_view, 'gt')
                compare_str = 'Show error';
            else
                compare_str = 'Show basis';
            end
            compare_enable = 'on';
        else
            compare_str = 'GT unavailable';
            compare_enable = 'off';
        end
        uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.84 0.955 0.14 0.036], 'String', compare_str, ...
            'BackgroundColor',[0.25 0.25 0.25], 'ForegroundColor','w', ...
            'FontSize',9, 'Enable', compare_enable, ...
            'Callback', @(~,~) toggleCompare());

        % ---- Shared slider (restored to current_slice) ----
        s = uicontrol('Style','slider', 'Units','normalized', ...
            'Position',[0.03 0.02 0.88 0.03], ...
            'Min',1,'Max',n_slices,'Value',current_slice, ...
            'SliderStep',[1/max(n_slices-1,1), 10/max(n_slices-1,1)], ...
            'BackgroundColor',[0.3 0.3 0.3]);

        lbl = uicontrol('Style','text', 'Units','normalized', ...
            'Position',[0.46 0.025 0.06 0.025], ...
            'String', sprintf('Slice %d / %d', current_slice, n_slices), ...
            'BackgroundColor','k', 'ForegroundColor','w', 'FontSize',10);

        s.Callback = @(~,~) updateSlice();
    end

    % ------------------------------------------------------------------
    % toggleMode  –  saves slice, flips mode, rebuilds layout, re-renders.
    % ------------------------------------------------------------------
    function toggleMode()
        current_slice = round(s.Value);
        if strcmpi(mode, 'combined')
            mode = 'separate';
        else
            mode = 'combined';
        end
        buildLayout();
        updateSlice();
    end

    function toggleCompare()
        if ~has_gt
            return
        end
        current_slice = round(s.Value);
        switch compare_view
            case 'basis'
                compare_view = 'gt';
            case 'gt'
                compare_view = 'error';
            otherwise
                compare_view = 'basis';
        end
        buildLayout();
        updateSlice();
    end

    function label = getModeLabel()
        switch compare_view
            case 'gt'
                label = sprintf('%s + GT compare', mode);
            case 'error'
                label = sprintf('%s + error maps', mode);
            otherwise
                label = mode;
        end
    end

    function makeCB(pos, cmap, vmin, vmax, label)
        ax = axes('Position', pos, 'Color','k', 'XColor','w', 'YColor','w');
        imagesc(ax, linspace(vmin, vmax, ncolors)');
        colormap(ax, cmap);
        ax.XTick = [];
        ax.YDir  = 'normal';
        ax.YTick = linspace(1, ncolors, 5);
        ax.YTickLabel = arrayfun(@(v) sprintf('%.6g', v), ...
            linspace(vmin, vmax, 5), 'UniformOutput', false);
        ax.TickDir = 'out';
        ax.YColor  = 'w';
        ylabel(ax, label, 'Color','w');
        ax.FontSize = 8;
    end

    function rgb = overlayPred(i_sl, gd_sl)
        i_n  = max(min((i_sl  - i_min)  / (i_max  - i_min),  1), 0);
        gd_n = max(min((gd_sl - gd_min) / (gd_max - gd_min), 1), 0);
        rgb = zeros([size(i_sl), 3]);
        rgb(:,:,3) = i_n;    % blue  = Iodine
        rgb(:,:,1) = gd_n;   % red   = Gadolinium
    end

    function rgb = overlayBasis(ib_sl, wb_sl)
        ib_n = max(min((ib_sl - ib_min) / (ib_max - ib_min), 1), 0);
        wb_n = max(min((wb_sl - wb_min) / (wb_max - wb_min), 1), 0);
        rgb = zeros([size(ib_sl), 3]);
        rgb(:,:,3) = ib_n;   % blue  = Iodine basis
        rgb(:,:,2) = wb_n;   % green = Water basis
    end

    function [ip, gp, ib, wb, igt, ggt, ierr, gerr, ttl] = getSlices(sl)
        switch lower(plane)
            case 'axial'
                ip = squeeze(i_vol(sl,:,:));   gp = squeeze(gd_vol(sl,:,:));
                ib = squeeze(ib_vol(sl,:,:));  wb = squeeze(wb_vol(sl,:,:));
                if has_gt
                    igt = squeeze(i_gt_vol(sl,:,:));
                    ggt = squeeze(gd_gt_vol(sl,:,:));
                    ierr = squeeze(i_err_vol(sl,:,:));
                    gerr = squeeze(gd_err_vol(sl,:,:));
                else
                    igt = [];
                    ggt = [];
                    ierr = [];
                    gerr = [];
                end
                ttl = sprintf('Axial  z = %d', sl);
            case 'coronal'
                ip = squeeze(i_vol(:,sl,:))';   gp = squeeze(gd_vol(:,sl,:))';
                ib = squeeze(ib_vol(:,sl,:))';  wb = squeeze(wb_vol(:,sl,:))';
                if has_gt
                    igt = squeeze(i_gt_vol(:,sl,:))';
                    ggt = squeeze(gd_gt_vol(:,sl,:))';
                    ierr = squeeze(i_err_vol(:,sl,:))';
                    gerr = squeeze(gd_err_vol(:,sl,:))';
                else
                    igt = [];
                    ggt = [];
                    ierr = [];
                    gerr = [];
                end
                ttl = sprintf('Coronal  x = %d', sl);
            case 'sagittal'
                ip = i_vol(:,:,sl)';   gp = gd_vol(:,:,sl)';
                ib = ib_vol(:,:,sl)';  wb = wb_vol(:,:,sl)';
                if has_gt
                    igt = i_gt_vol(:,:,sl)';
                    ggt = gd_gt_vol(:,:,sl)';
                    ierr = i_err_vol(:,:,sl)';
                    gerr = gd_err_vol(:,:,sl)';
                else
                    igt = [];
                    ggt = [];
                    ierr = [];
                    gerr = [];
                end
                ttl = sprintf('Sagittal  y = %d', sl);
        end

    ip = rot90(fliplr(ip));
    gp = rot90(fliplr(gp));
    ib = rot90(fliplr(ib));
    wb = rot90(fliplr(wb));
    if has_gt
        igt = rot90(fliplr(igt));
        ggt = rot90(fliplr(ggt));
        ierr = rot90(fliplr(ierr));
        gerr = rot90(fliplr(gerr));
    end

    end

    function rgb = scalarToRGB(data, vmin, vmax, cmap)
        data_n = max(0, min(1, (data - vmin) / (vmax - vmin)));
        idx    = gray2ind(data_n, ncolors);
        rgb    = ind2rgb(idx, cmap);
    end

    function updateSlice()
        sl = round(s.Value);
        lbl.String = sprintf('Slice %d / %d', sl, n_slices);
        [ip, gp, ib, wb, igt, ggt, ierr, gerr, ttl] = getSlices(sl);

        if strcmpi(mode, 'combined')
            imshow(overlayPred(ip,  gp),  'Parent', ax_pred);
            if strcmp(compare_view, 'gt')
                imshow(overlayPred(igt, ggt), 'Parent', ax_basis);
            elseif strcmp(compare_view, 'error')
                imshow(overlayPred(ierr, gerr), 'Parent', ax_basis);
            else
                imshow(overlayBasis(ib, wb),  'Parent', ax_basis);
            end
            title(ax_pred,  ['Predicted CAs (I/Gd)–  ' ttl], 'Color','w','FontSize',11);
            if strcmp(compare_view, 'gt')
                title(ax_basis, ['Ground truth CAs (I/Gd) – ' ttl], 'Color','w','FontSize',11);
            elseif strcmp(compare_view, 'error')
                title(ax_basis, ['Absolute error |prediction - GT| (I/Gd) – ' ttl], 'Color','w','FontSize',11);
            else
                title(ax_basis, ['Basis pair (I/Wa) – ' ttl], 'Color','w','FontSize',11);
            end
        else
            imshow(scalarToRGB(ip, i_min,  i_max,  cmap_blue),  'Parent', ax_ip);
            imshow(scalarToRGB(gp, gd_min, gd_max, cmap_red),   'Parent', ax_gp);
            if strcmp(compare_view, 'gt')
                imshow(scalarToRGB(igt, i_min,  i_max,  cmap_blue), 'Parent', ax_ib);
                imshow(scalarToRGB(ggt, gd_min, gd_max, cmap_red),  'Parent', ax_wb);
            elseif strcmp(compare_view, 'error')
                imshow(scalarToRGB(ierr, 0, i_max, cmap_blue), 'Parent', ax_ib);
                imshow(scalarToRGB(gerr, 0, gd_max, cmap_red), 'Parent', ax_wb);
            else
                imshow(scalarToRGB(ib, ib_min, ib_max, cmap_blue),  'Parent', ax_ib);
                imshow(scalarToRGB(wb, wb_min, wb_max, cmap_green), 'Parent', ax_wb);
            end
            title(ax_ip, ['Predicted Iodine – '     ttl], 'Color','w','FontSize',10);
            title(ax_gp, ['Predicted Gadolinium – ' ttl], 'Color','w','FontSize',10);
            if strcmp(compare_view, 'gt')
                title(ax_ib, ['Ground truth Iodine – '     ttl], 'Color','w','FontSize',10);
                title(ax_wb, ['Ground truth Gadolinium – ' ttl], 'Color','w','FontSize',10);
            elseif strcmp(compare_view, 'error')
                title(ax_ib, ['Absolute iodine error – ' ttl], 'Color','w','FontSize',10);
                title(ax_wb, ['Absolute gadolinium error – ' ttl], 'Color','w','FontSize',10);
            else
                title(ax_ib, ['Basis Iodine – ' ttl], 'Color','w','FontSize',10);
                title(ax_wb, ['Basis Water – '  ttl], 'Color','w','FontSize',10);
            end
        end
    end

end