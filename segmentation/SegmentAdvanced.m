function [PE0, PVC0, DeltaPE, DeltaPVC, Mono70keV, MaskTot] = Segment(PE_ref, PVC_ref, PE, PVC, ...,
    PEraw, PVCraw, segmentation, labels, organs_to_segment, ConstantOnOrgans)

%% Params

kernelsize = 5;

%% Get monos

AirThreshold = -50;
diffMaxHUThreshold = 250;
E = 40;
MuPE  = 0.1*GetMu('polyethylene', E);
MuPVC = 0.1*GetMu('pvc', E);

Mono70keV = PE * MuPE + PVC * MuPVC; 
Mono70keVHU = TransformHU(Mono70keV, E);

Mono70keV_VNC = PE_ref * MuPE + PVC_ref * MuPVC; 
Mono70keVHU_VNC = TransformHU(Mono70keV_VNC, E);

diffMask = imgaussfilt(Mono70keVHU - Mono70keVHU_VNC, 8) / diffMaxHUThreshold;
diffMask(diffMask > 1) = 1;
diffMask(diffMask< 0) = 0;

AirMask = imgaussfilt(Mono70keVHU, 2) > AirThreshold;

%% Total mask (where we expect non-zero values)

MaskTot = zeros(size(PE));

%% Labels

Norgans = numel(labels);

%% Get PE/PVC maps

PE0 = PE;
PVC0 = PVC;

dPE = PE - PE_ref;
dPVC = PVC - PVC_ref;

DeltaPE = zeros(size(PE));
DeltaPVC = zeros(size(PVC));

if any(strcmp(organs_to_segment, 'all'))
    DeltaPE = PE - PE_ref;
    DeltaPVC = PVC - PVC_ref;
else
    for i=1:Norgans
        organ = labels{i};

        if any(strcmp(organs_to_segment, organ))

           mask = zeros(size(segmentation));
           mask(segmentation==i) = 1;
           mask = mask .* AirMask .* diffMask;

           MaskTot = MaskTot + mask;

           %% Iterate over submasks

           CC = bwconncomp(mask);
           numRegions = CC.NumObjects;

           for k = 1:numRegions
               submask = false(size(mask));
               submask(CC.PixelIdxList{k}) = true;

               if sum(submask, 'all') > 10

                   organAvgPE = mean(nonzeros(submask.*PE_ref), 'all');
                   organAvgPVC = mean(nonzeros(submask.*PVC_ref), 'all');


                   maskPEBackground = submask * organAvgPE;
                   maskPVCBackground = submask * organAvgPVC;

                   maskPE = PE .*submask;
                   maskPVC = PVC.*submask;

                   maskPE0 = PE0 .*submask;
                   maskPVC0 = PVC0 .* submask;

                   maskPEraw = PEraw .* submask;
                   maskPVCraw = PVCraw .*submask;

                   if ConstantOnOrgans
                       maskPE = mean(nonzeros(maskPE), 'all') * submask;
                       maskPVC = mean(nonzeros(maskPVC), 'all') * submask;
                   else
                       maskPEBackground = maskedLowPass(maskPEBackground, submask, kernelsize);
                       maskPVCBackground = maskedLowPass(maskPVCBackground, submask, kernelsize);
                       maskPE = maskedLowPass(maskPE, submask, kernelsize);
                       maskPVC = maskedLowPass(maskPVC, submask, kernelsize);
                       maskPE0 = maskedLowPass(maskPE0, submask, kernelsize);
                       maskPVC0 = maskedLowPass(maskPVC0, submask, kernelsize);
                   end



                   maskDeltaPE = maskPE - maskPEBackground;
                   maskDeltaPVC = maskPVC - maskPVCBackground;
                   maskDeltaPEraw = maskPEraw - maskPEBackground;
                   maskDeltaPVCraw = maskPVCraw - maskPVCBackground;

                   %% Ensure Delta PE / Delta PVC are not averaged out

                   maskDeltaPE = maskDeltaPE * mean(nonzeros(maskDeltaPEraw), 'all') / mean(nonzeros(maskDeltaPE), 'all');
                   maskDeltaPVC = maskDeltaPVC * mean(nonzeros(maskDeltaPVCraw), 'all') / mean(nonzeros(maskDeltaPVC), 'all');

                   % Updata PE/PVC maps

                   DeltaPE = DeltaPE + maskDeltaPE;
                   DeltaPVC = DeltaPVC + maskDeltaPVC;

                   PE0(submask) = 0;
                   PVC0(submask) = 0;
                   PE0 = PE0 + maskPE0 - maskDeltaPE;
                   PVC0 = PVC0 + maskPVC0 - maskDeltaPVC;
               end
           end
           end
        end
    end

% Normalize MaskTot

MaskTot = MaskTot > 0;
end


function M = maskedLowPass(Mk, mask, kernelsize)
    h = ones(kernelsize, kernelsize);
    num = conv2(Mk, h, 'same');
    denom = conv2(mask, h, 'same');
    denom(denom==0) = NaN;
    M = num./denom;
    M(~mask) = 0;

end

function MonoHU = TransformHU(Mono, E)
    MuWa = GetMu('water', E);
    MuAir = GetMu('air', E);
    MonoHU = 1000 * (Mono - 0.1*MuWa) / (0.1*MuWa - 0.1*MuAir);
end