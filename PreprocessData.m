function [AB, AR, AC, MaskTot, detector_cfg] = PreprocessData(a, aRef, seg, labels, pixelstomm, phantom_type, visualization)


if nargin < 6
   visualization = 0; 
end

%% Apply low pas filter

PE =  a{1};
PVC = a{2};



options = struct( ...
    'SearchWindowSize', 11, ...     % Increase to look farther for similar patches
    'ComparisonWindowSize', 11, ... % Larger patch to match texture
    'DegreeOfSmoothing', 8);        % Increase to smooth more aggressively

PE_filt = imnlmfilt(PE, options);
PVC_filt = imnlmfilt(PVC, options);



%% segmentation

switch phantom_type
    case 'antropomorphic'
        
        PE_ref =  aRef{1};
        PVC_ref = aRef{2};
        PE_ref_filt = imnlmfilt(PE_ref, options);
        PVC_ref_filt = imnlmfilt(PVC_ref, options);
        
        organs_to_segment = {'stomach', 'liver', 'kidney_left', 'kidney_right', ...
            'stomach', 'aorta', 'inferior_vena_cava', 'colon', 'small_bowel', 'spleen'};
        % organs_to_segment = {'aorta', 'inferior_vena_cava'};
        % organs_to_segment = {'liver', 'kidney_left', 'kidney_right'};

        organs_to_segment = {'small_bowel', 'blood', 'kidney_left', ...
           'kidney_right', 'pancreas', 'inferior_vena_cava', 'aorta', ..., 
           'stomach', 'adrenal_gland_right', 'adrenal_gland_left', ...,
           'duodenum', 'liver', 'portal_vein_and_splenic_vein', 'colon'};

        % organs_to_segment = {'liver'};

        ConstantOnOrgans = 0;
        % 
        [PE, PVC, DeltaPE, DeltaPVC, ~, MaskTot] =  SegmentAdvanced(PE_ref_filt, PVC_ref_filt, ...,
            PE_filt, PVC_filt, PE, PVC ...  
            , seg, [labels], organs_to_segment, ConstantOnOrgans);
   
        
    case 'gammex'
        I = a{3};
        Gd = a{4};
        [PE, PVC, DeltaPE, DeltaPVC, ~, MaskTot] = SegmentSimple(PE_filt, PVC_filt, I, Gd);
end

%% Optional: Visualize Segmentations

if visualization
    plot_input({PE, PVC}, {DeltaPE, DeltaPVC});
end


%% Calculate sinograms

% GT sinogtrams

%% Fanbeam

load('bowtiecorrection_gammas.mat', 'Detector_Gammas')

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

    Idx = abs(gamma) < 0.8 * max(abs(Detector_Gammas));
    gamma = gamma(Idx);
    
    s = detector_cfg.pixelstomm * s;
    s = s(Idx, :);
    A = [A, s];
end


% Segmented sinograms

PE_sino = pixelstomm*fanbeam(PE, D_px, ...,
    'FanSensorGeometry', 'arc', ...
    'FanSensorSpacing', spacing, ...
    'FanRotationIncrement', rot_increment);

PVC_sino = pixelstomm*fanbeam(PVC, D_px, ...,
    'FanSensorGeometry', 'arc', ...
    'FanSensorSpacing', spacing, ...
    'FanRotationIncrement', rot_increment);

DeltaPE_sino = pixelstomm*fanbeam(DeltaPE, D_px, ...,
    'FanSensorGeometry', 'arc', ...
    'FanSensorSpacing', spacing, ...
    'FanRotationIncrement', rot_increment);
DeltaPVC_sino = pixelstomm*fanbeam(DeltaPVC, D_px,  ...,
    'FanSensorGeometry', 'arc', ...
    'FanSensorSpacing', spacing, ...
    'FanRotationIncrement', rot_increment);

PE_sino = PE_sino(Idx, :);
PVC_sino = PVC_sino(Idx, :);
DeltaPE_sino = DeltaPE_sino(Idx, :);
DeltaPVC_sino = DeltaPVC_sino(Idx, :);

%% Add bowtie Correction

load('bowtiecorrection.mat', 'Correction');

CorrectionPE = Correction(:, 1);
CorrectionPVC = Correction(:, 2);

idx_est = interp1(Detector_Gammas, 1:length(Detector_Gammas), gamma, 'nearest');

CorrectionPE = CorrectionPE(idx_est);
CorrectionPVC = CorrectionPVC(idx_est);
CorrectionGammas = Detector_Gammas(idx_est);

M = size(PE_sino, 2);
CorrectionPE_sino = repmat(CorrectionPE, 1, M);
CorrectionPVC_sino = repmat(CorrectionPVC, 1, M);

PE_sino = PE_sino + CorrectionPE_sino;
PVC_sino = PVC_sino + CorrectionPVC_sino;




AB = {PE_sino, PVC_sino};                   % Material Basis sinograms without contrast
AR = {DeltaPE_sino, DeltaPVC_sino};         % Contrast Contribution to basis sinograms
AC = {A{3}, A{4}};                          % GT Contrast sinograms

%% Optional: Plot segmentation sinograms

if visualization
   plot_sinos(AB, AR, rotAngles, gamma)
end


end