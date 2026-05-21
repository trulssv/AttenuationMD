function [PE0, PVC0, DeltaPE, DeltaPVC, mono70keV, MaskTot] = SegmentSimple(PE, PVC, I , Gd)

%% Segmenation


MaskTot = max(abs(I), 3 * abs(Gd));
MaskTot(MaskTot <5e-5) = 0;
MaskTot(MaskTot ~= 0) = 1;

%% Air mask
E = 75;
muPE = 0.1 * GetMu('polyethylene', E);
muPVC = 0.1 * GetMu('pvc', E);
mono70keV = PE * muPE + PVC * muPVC;
mono70keV = TransformHU(mono70keV, E);
mono70keV = imgaussfilt(mono70keV, 2);
airmask = mono70keV > -100;

MaskTot = MaskTot .* airmask;

%% Perform segmentation

meanPE = mean(nonzeros(airmask .* ~MaskTot .* PE), 'all');
meanPVC = mean(nonzeros(airmask .* ~MaskTot .* PVC), 'all');
PEvnc = ~MaskTot .* PE + MaskTot*meanPE;
PVCvnc = ~MaskTot .* PVC + MaskTot*meanPVC;

CC = bwconncomp(MaskTot);
numRegions = CC.NumObjects;

DeltaPE = zeros(size(MaskTot));
DeltaPVC = zeros(size(MaskTot));

for k = 1:numRegions
    submask = false(size(MaskTot));
    submask(CC.PixelIdxList{k}) = true;

    if sum(submask, 'all') > 100

        DeltaPEsubmask = submask .*(PE- PEvnc);
        DeltaPVCsubmask = submask .*(PVC- PVCvnc);

        DeltaPEsubmask_filt = maskedLowPass(DeltaPEsubmask, submask, 3);
        DeltaPVCsubmask_filt = maskedLowPass(DeltaPVCsubmask, submask, 3);

        DeltaPE = DeltaPE +DeltaPEsubmask_filt;
        DeltaPVC = DeltaPVC + DeltaPVCsubmask_filt;
    else
        MaskTot(submask) = 0;
    end
end


%% Apply low pas filter
options = struct( ...
    'SearchWindowSize', 15, ...      % Increase to look farther for similar patches
    'ComparisonWindowSize', 15, ...   % Larger patch to match texture
    'DegreeOfSmoothing', 10);        % Increase to smooth more aggressively

PE0 = PE;
PVC0 = PVC;

PE0 = imnlmfilt(PE0, options);
PVC0 = imnlmfilt(PVC0, options);