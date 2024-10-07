function [] = ProcessRawDataFiles_FP_EEG(rawDataFiles)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the force sensor and neural bands. Create a threshold for binarized movement/whisking if
%            one does not already exist.
%________________________________________________________________________________________________________________________

% Raw data file analysis
for a = 1:size(rawDataFiles,1)
    rawDataFile = rawDataFiles(a,:);
    disp(['Analyzing RawData file ' num2str(a) ' of ' num2str(size(rawDataFiles, 1)) '...']); disp(' ')
    [animalID, fileDate, fileID] = GetFileInfo_FP(rawDataFile);
    strDay = ConvertDate_FP(fileDate);
    procDataFile = ([animalID '_' fileID '_ProcData.mat']);
    disp(['Generating ' procDataFile '...']); disp(' ')
    load(rawDataFile);
    % Transfer RawData notes to ProcData structure.
    ProcData.notes = RawData.notes;
    ProcData.notes.dsFs = 30; % downsampled Fs
    analogExpectedLength = ProcData.notes.trialDuration_sec*ProcData.notes.analogSamplingRate;
    dsExpectedLength = ProcData.notes.trialDuration_sec*ProcData.notes.dsFs;
    %% identify the solenoids by amplitude.
    ProcData.data.stimulations.LPadSol = find(diff(RawData.data.stimulations) == 1)/RawData.notes.analogSamplingRate;
    ProcData.data.stimulations.RPadSol = find(diff(RawData.data.stimulations) == 2)/RawData.notes.analogSamplingRate;
    ProcData.data.stimulations.AudSol = find(diff(RawData.data.stimulations) == 3)/RawData.notes.analogSamplingRate;
    %% resample doric data to 30 Hz
%     ProcData.data.Rhodamine.LH = (resample(RawData.data.LH.dF.F560,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
%     ProcData.data.Rhodamine.LH = ProcData.data.Rhodamine.LH(1:dsExpectedLength);
    ProcData.data.Rhodamine.RH = (resample(RawData.data.RH.dF.F560,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.Rhodamine.RH = ProcData.data.Rhodamine.RH(1:dsExpectedLength);
%     ProcData.data.GFP.LH = (resample(RawData.data.LH.dF.F465,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
%     ProcData.data.GFP.LH = ProcData.data.GFP.LH(1:dsExpectedLength);
    ProcData.data.GFP.RH = (resample(RawData.data.RH.dF.F465,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.GFP.RH = ProcData.data.GFP.RH(1:dsExpectedLength);
    %% process neural data into its various forms.
    neuralDataTypes = {'EEG_LH','EEG_RH'};
    for c = 1:length(neuralDataTypes)
        neuralDataType = neuralDataTypes{1,c};
        % process power of the frequency and the raw signal
        % EEG Band [1 - 100]. Captures the total EEG frequency spectrum
        [EEGPower,EEGSignal,~] = ProcessNeuro_FP_EEG_Signal(RawData,analogExpectedLength,'EEG',neuralDataType);
        ProcData.data.(neuralDataType).EEGPower = EEGPower;
        ProcData.data.(neuralDataType).EEGSignal = EEGSignal;
        % Gamma Band [40 - 100]
        [gammaBandPower,gammaBandSignal,~] = ProcessNeuro_FP_EEG_Signal(RawData,analogExpectedLength,'Gam',neuralDataType);
        ProcData.data.(neuralDataType).gammaBandPower = gammaBandPower;
        ProcData.data.(neuralDataType).gammaBandSignal = gammaBandSignal;
        % Beta [13 - 30 Hz]
        [betaBandPower,betaBandSignal,~] = ProcessNeuro_FP_EEG_Signal(RawData,analogExpectedLength,'Beta',neuralDataType);
        ProcData.data.(neuralDataType).betaBandPower = betaBandPower;
        ProcData.data.(neuralDataType).betaBandSignal = betaBandSignal;
        % Alpha [8 - 12 Hz]
        [alphaBandPower,alphaBandSignal,~] = ProcessNeuro_FP_EEG_Signal(RawData,analogExpectedLength,'Alpha',neuralDataType);
        ProcData.data.(neuralDataType).alphaBandPower = alphaBandPower;
        ProcData.data.(neuralDataType).alphaBandSignal = alphaBandSignal;
        % Theta [4 - 8 Hz]
        [thetaBandPower,thetaBandSignal,~] = ProcessNeuro_FP_EEG_Signal(RawData,analogExpectedLength,'Theta',neuralDataType);
        ProcData.data.(neuralDataType).thetaBandPower = thetaBandPower;
        ProcData.data.(neuralDataType).thetaBandSignal = thetaBandSignal;
        % Delta [1 - 4 Hz]
        [deltaBandPower,deltaBandSignal,~] = ProcessNeuro_FP_EEG_Signal(RawData,analogExpectedLength,'Delta',neuralDataType);
        ProcData.data.(neuralDataType).deltaBandPower = deltaBandPower;
        ProcData.data.(neuralDataType).deltaBandSignal = deltaBandSignal;
    end
    %% Patch and binarize the whisker angle and set the resting angle to zero degrees.
    filtThreshold = 20;
    filtOrder = 2;
    [z,p,k] = butter(filtOrder,filtThreshold/(RawData.notes.whiskCamSamplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    filteredWhiskers = filtfilt(sos,g,RawData.data.whiskerAngle - mean(RawData.data.whiskerAngle));
    resampledWhisk = resample(filteredWhiskers,ProcData.notes.dsFs,RawData.notes.whiskCamSamplingRate);
    % Binarize the whisker waveform (wwf)
    threshfile = dir('*_Thresholds.mat');
    if ~isempty(threshfile)
        load(threshfile.name)
    end
    [ok] = CheckForThreshold_FP(['binarizedWhiskersLower_' strDay],animalID);
    if ok == 0
        [whiskersThresh1,whiskersThresh2] = CreateWhiskThreshold_FP(resampledWhisk,ProcData.notes.dsFs);
        Thresholds.(['binarizedWhiskersLower_' strDay]) = whiskersThresh1;
        Thresholds.(['binarizedWhiskersUpper_' strDay]) = whiskersThresh2;
        save([animalID '_Thresholds.mat'], 'Thresholds');
    end
    load([animalID '_Thresholds.mat']);
    binWhisk = BinarizeWhiskers_FP(resampledWhisk,ProcData.notes.dsFs,Thresholds.(['binarizedWhiskersLower_' strDay]),Thresholds.(['binarizedWhiskersUpper_' strDay]));
    [linkedBinarizedWhiskers] = LinkBinaryEvents_FP(gt(binWhisk,0),[round(ProcData.notes.dsFs/3),0]);
    inds = linkedBinarizedWhiskers == 0;
    restAngle = mean(resampledWhisk(inds));
    ProcData.data.whiskerAngle = resampledWhisk - restAngle;
    ProcData.data.binWhiskerAngle = binWhisk;
    %% Downsample and binarize the force sensor.
    trimmedForce = RawData.data.forceSensor(1:min(analogExpectedLength,length(RawData.data.forceSensor)));
    % Filter then downsample the Force Sensor waveform to desired frequency
    filtThreshold = 20;
    filtOrder = 2;
    [z,p,k] = butter(filtOrder,filtThreshold/(ProcData.notes.analogSamplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    filtForceSensor = filtfilt(sos,g,trimmedForce);
    ProcData.data.forceSensor = resample(filtForceSensor,ProcData.notes.dsFs,ProcData.notes.analogSamplingRate);
    % Binarize the force sensor waveform
    threshfile = dir('*_Thresholds.mat');
    if ~isempty(threshfile)
        load(threshfile.name)
    end
    [ok] = CheckForThreshold_FP(['binarizedForceSensor_' strDay],animalID);
    if ok==0
        [forceSensorThreshold] = CreateForceSensorThreshold_FP_Up(ProcData.data.forceSensor);
        Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
        save([animalID '_Thresholds.mat'],'Thresholds');
    end
    if ok==1
            TR = load([animalID '_Thresholds.mat']);
            if TR.Thresholds.(['binarizedForceSensor_' strDay])==0.005 || TR.Thresholds.(['binarizedForceSensor_' strDay])==0.03
                    [forceSensorThreshold] = CreateForceSensorThreshold_FP_Up(ProcData.data.forceSensor);
                    Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
                    save([animalID '_Thresholds.mat'],'Thresholds');
            end
    end

    ProcData.data.binForceSensor = BinarizeForceSensor_FP(ProcData.data.forceSensor,Thresholds.(['binarizedForceSensor_' strDay]));
    %% EMG    
    trimmedEMG = RawData.data.EMG(1:min(analogExpectedLength,length(RawData.data.EMG)));

    fpass = [300,3000];
    [z,p,k] = butter(3,fpass/(ProcData.notes.analogSamplingRate/2));
    [sos,g] = zp2sos(z,p,k);
    filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
    kernelWidth = 0.5;
    smoothingKernel = gausswin(kernelWidth*ProcData.notes.analogSamplingRate)/sum(gausswin(kernelWidth*ProcData.notes.analogSamplingRate));
    EMGPwr = log10(conv(filtEMG.^2,smoothingKernel,'same'));
    resampEMG = resample(EMGPwr,ProcData.notes.dsFs,ProcData.notes.analogSamplingRate);
    ProcData.data.EMG.power = resampEMG;

    fpass = [10,100];
    [z,p,k] = butter(3,fpass/(ProcData.notes.analogSamplingRate/2));
    [sos,g] = zp2sos(z,p,k);
    filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
    resampEMG = resample(filtEMG,ProcData.notes.dsFs,ProcData.notes.analogSamplingRate);
    ProcData.data.EMG.emg = resampEMG;
    %% Pupil
    ProcData.data.Pupil.pupilArea = RawData.data.pupilArea; 
    ProcData.data.Pupil.mmarea = RawData.data.pupilmmArea;
    ProcData.data.Pupil.Diameter = RawData.data.pupilDiameter;
    ProcData.data.Pupil.mmDiameter = RawData.data.pupilmmDiameter;

    ProcData.data.Pupil.Major = RawData.data.pupilpatchMajor;
    ProcData.data.Pupil.Minor = RawData.data.pupilpatchMinor;
    ProcData.data.Pupil.CentroidX = RawData.data.pupilpatchCentroidX;
    ProcData.data.Pupil.CentroidY = RawData.data.pupilpatchCentroidY;

   
    ProcData.data.Pupil.blinkFrames = RawData.data.Pupil.blinkFrames;
    ProcData.data.Pupil.blinkInds = RawData.data.Pupil.blinkInds;
    ProcData.data.Pupil.eyeROI = RawData.data.Pupil.eyeROI;
    ProcData.data.Pupil.firstFrame = RawData.data.Pupil.firstFrame;
    %% save the processed data
    save(procDataFile,'ProcData')
end

end
