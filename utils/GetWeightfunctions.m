function Weights = GetWeightfunctions(eActVectorkeV)

addpath /media/data/trulssv/IodineGadolinium/IodineGadolinium/CatSim/Git_20240312/main_GE/base/materials/
addpath /media/data/trulssv/IodineGadolinium/IodineGadolinium/CatSim/Git_20240312/main_GE/base/mfiles/
addpath /media/data/trulssv/IodineGadolinium/IodineGadolinium/CatSim/2024_06_13_MCM_BMD_code
simulatedResponseFunctionPathAndFilename = '/media/data/trulssv/IodineGadolinium/IodineGadolinium/CatSim/pluto_response_database/response_matrix/Pluto_spectral_response_withChargeSharingGaussian_pixel_50.mat';
% simulatedResponseFunctionPathAndFilename = '/media/data/trulssv/IodineGadolinium/IodineGadolinium/CatSim/pluto_response_database/response_matrix/Pluto_spectral_response_withChargeSharingUniform_pixel_50.mat';


kVp = 120;


eNoiseSigmakeV = 2.75; %Measurement by Bj?rn Cederstr?m

prepatientPhotonsPerDel= 1e6;

%Tube spectrum, filtering
additionalTubeFiltermmAl=0; %3.25 mm additional filtering gives HVL=7.4 mm Al, which agrees with the large bowtie. 

nSkippedEActStepsInBeginningOfRespFunc=0; %Don't skip any low energies for now

[originalTubeSpectrum,originalTubeSpectrumEnergies] = Spectrum_Read(sprintf('/media/data/trulssv/IodineGadolinium/IodineGadolinium/Git_20250113/PCCT/spectrum/xspect_tar10.5_%d_filtNom.dat',kVp));
unfilteredTubeSpectrum=interp1(originalTubeSpectrumEnergies,originalTubeSpectrum,eActVectorkeV,'linear',0);

muSimm_1=0.1*GetMu('si',eActVectorkeV);
muAlmm_1=0.1*GetMu('al',eActVectorkeV);

filteredTubeSpectrumNormalizedToOne = unfilteredTubeSpectrum.*exp(-muAlmm_1*additionalTubeFiltermmAl); %Changes on the next line!
filteredTubeSpectrumNormalizedToOne = filteredTubeSpectrumNormalizedToOne/sum(filteredTubeSpectrumNormalizedToOne);


deadSiLayermm=0.285; %TODO check if this is up to date
depthSegmentLengthsmm = [9.2 18.8 9.0]; %Should apply to one of the central modules. TODO check if this is up to date
nDepthSegments=length(depthSegmentLengthsmm);
detectionEfficiencyInDepthSegments = zeros(nDepthSegments,length(eActVectorkeV));
for depthSegmentNo=1:nDepthSegments
    detectionEfficiencyInDepthSegments(depthSegmentNo,:) = exp(-muSimm_1*(deadSiLayermm+sum(depthSegmentLengthsmm(1:depthSegmentNo-1)))).*(1-exp(-muSimm_1*depthSegmentLengthsmm(depthSegmentNo)));
end

simulatedEnergyResponse = load(simulatedResponseFunctionPathAndFilename,'Evec0','Dvec0','res_mat');
energyResponseFunctionForDepthSegments = repmat(permute(simulatedEnergyResponse.res_mat(nSkippedEActStepsInBeginningOfRespFunc+1:end,:),[2 3 1]),[1 nDepthSegments 1]); %Use the same response function for each depth segment for now
eActVectorForSimulatedEnergyResponse=simulatedEnergyResponse.Evec0(nSkippedEActStepsInBeginningOfRespFunc+1:end);
binThresholdsForDepthSegmentskeV = repmat(15:10:85,[nDepthSegments 1]); %Dummy, TODO replace by real thresholds once we figure out what they are.
depthSegmentsToUse=ones(1,nDepthSegments); %Possible to skip some depth segments for example if they are not working
aggregateDepthSegments=false;


filteredTubeSpectrum=prepatientPhotonsPerDel*filteredTubeSpectrumNormalizedToOne;
Weights = getForwardModelWeightFunction(filteredTubeSpectrum, detectionEfficiencyInDepthSegments, eActVectorkeV, energyResponseFunctionForDepthSegments, eActVectorForSimulatedEnergyResponse, simulatedEnergyResponse.Dvec0, binThresholdsForDepthSegmentskeV, depthSegmentsToUse, eNoiseSigmakeV,aggregateDepthSegments);