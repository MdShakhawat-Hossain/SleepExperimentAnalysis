function [] = ProcessRawDataFiles_FP_GRABNE(rawDataFiles)
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
    if exist(procDataFile,"file") == 2
        load(procDataFile);
    end

    % Transfer RawData notes to ProcData structure.
    ProcData.notes = RawData.notes;
    ProcData.notes.dsFs = 30; % downsampled Fs
    analogExpectedLength = ProcData.notes.trialDuration_sec*ProcData.notes.analogSamplingRate;
    dsExpectedLength = ProcData.notes.trialDuration_sec*ProcData.notes.dsFs;
    %% identify the solenoids by amplitude.
    ProcData.data.stimulations.LPadSol = find(diff(RawData.data.stimulations) == 1)/RawData.notes.analogSamplingRate;
    ProcData.data.stimulations.RPadSol = find(diff(RawData.data.stimulations) == 2)/RawData.notes.analogSamplingRate;
    ProcData.data.stimulations.AudSol = find(diff(RawData.data.stimulations) == 3)/RawData.notes.analogSamplingRate;
    ProcData.data.stimulations.OptoStim = find(diff(RawData.data.stimulations) == 4)/RawData.notes.analogSamplingRate;
    %% resample fiber data to 30 Hz
        % remove some previous data
    if isfield(ProcData.data,'GFP') == 1
    ProcData.data = rmfield(ProcData.data,'GFP');
    end
    if isfield(ProcData.data,'CBV') == 1
    ProcData.data = rmfield(ProcData.data,'CBV');
    end
    if isfield(ProcData.data,'CBV') == 1
    ProcData.data = rmfield(ProcData.data,'CBV');
    end
    if isfield(ProcData.data,'Isos') == 1
    ProcData.data = rmfield(ProcData.data,'Isos');
    end

    % ZScored
    ProcData.data.CBV.Z_ACh = (resample(RawData.data.ACh.dFF0_z.F560,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.CBV.Z_ACh = ProcData.data.CBV.Z_ACh(1:dsExpectedLength);
    ProcData.data.CBV.Z_NE = (resample(RawData.data.NE.dFF0_z.F560,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.CBV.Z_NE = ProcData.data.CBV.Z_NE(1:dsExpectedLength);

    ProcData.data.GFP.Z_ACh = (resample(RawData.data.ACh.dFF0_z.F465,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.GFP.Z_ACh = ProcData.data.GFP.Z_ACh(1:dsExpectedLength);
    ProcData.data.GFP.Z_NE = (resample(RawData.data.NE.dFF0_z.F465,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.GFP.Z_NE = ProcData.data.GFP.Z_NE(1:dsExpectedLength);

    ProcData.data.Isos.Z_ACh = (resample(RawData.data.ACh.dFF0_z.F405,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.Isos.Z_ACh = ProcData.data.Isos.Z_ACh(1:dsExpectedLength);
    ProcData.data.Isos.Z_NE = (resample(RawData.data.NE.dFF0_z.F405,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.Isos.Z_NE = ProcData.data.Isos.Z_NE(1:dsExpectedLength);

    % Percentage
    ProcData.data.CBV.P_ACh = (resample(RawData.data.ACh.dFF0_p.F560,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.CBV.P_ACh = ProcData.data.CBV.P_ACh(1:dsExpectedLength);
    ProcData.data.CBV.P_NE = (resample(RawData.data.NE.dFF0_p.F560,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.CBV.P_NE = ProcData.data.CBV.P_NE(1:dsExpectedLength);

    ProcData.data.GFP.P_ACh = (resample(RawData.data.ACh.dFF0_p.F465,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.GFP.P_ACh = ProcData.data.GFP.P_ACh(1:dsExpectedLength);
    ProcData.data.GFP.P_NE = (resample(RawData.data.NE.dFF0_p.F465,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.GFP.P_NE = ProcData.data.GFP.P_NE(1:dsExpectedLength);

    ProcData.data.Isos.P_ACh = (resample(RawData.data.ACh.dFF0_p.F405,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.Isos.P_ACh = ProcData.data.Isos.P_ACh(1:dsExpectedLength);
    ProcData.data.Isos.P_NE = (resample(RawData.data.NE.dFF0_p.F405,ProcData.notes.dsFs,ProcData.notes.TDT.DataFs));
    ProcData.data.Isos.P_NE = ProcData.data.Isos.P_NE(1:dsExpectedLength);
    %% process neural data into its various forms.
    neuralDataTypes = {'cortical_LH','cortical_RH'};%,'cortical_RH','hippocampus'};
    for c = 1:length(neuralDataTypes)
        neuralDataType = neuralDataTypes{1,c};
        %median filter the data to remove sudden spikes
        RawData.data.(neuralDataType) = medfilt1(RawData.data.(neuralDataType),10); % 10th order median filter

        % 60HZ notch filter
        [zN, pN, kN] = butter(3, [59 61]./(RawData.notes.analogSamplingRate/2), 'stop'); % 3rd order band stop filter
        [sosN,gN]=zp2sos(zN,pN,kN); 
        RawData.data.(neuralDataType) = filtfilt(sosN,gN,RawData.data.(neuralDataType)); % notch filter the raw data
        RawData.data.(neuralDataType)(1:RawData.notes.analogSamplingRate) = RawData.data.(neuralDataType)(RawData.notes.analogSamplingRate+1:2*RawData.notes.analogSamplingRate); 
        % process power of the frequency and the raw signal
        % MUA Band [300 - 3000]
%         [muaPower,muaSignal,~] = ProcessNeuro_FP_EEG_Signal(RawData,analogExpectedLength,'MUA',neuralDataType);
%         ProcData.data.(neuralDataType).muaPower = muaPower;
%         ProcData.data.(neuralDataType).muaSignal = muaSignal;
        % cortical Band [1 - 100]. Captures the total cortical frequency spectrum
        [corticalPower,corticalSignal,~] = ProcessNeuro_FP_EEG_Signal(RawData,analogExpectedLength,'EEG',neuralDataType);
        ProcData.data.(neuralDataType).corticalPower = corticalPower;
        ProcData.data.(neuralDataType).corticalSignal = corticalSignal;
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
            if  TR.Thresholds.(['binarizedForceSensor_' strDay])==0.03 % TR.Thresholds.(['binarizedForceSensor_' strDay])==0.005 ||
                    [forceSensorThreshold] = CreateForceSensorThreshold_FP_Up(ProcData.data.forceSensor);
                    Thresholds.(['binarizedForceSensor_' strDay]) = forceSensorThreshold;
                    save([animalID '_Thresholds.mat'],'Thresholds');
            end
    end

    ProcData.data.binForceSensor = BinarizeForceSensor_FP(ProcData.data.forceSensor,Thresholds.(['binarizedForceSensor_' strDay]));
    %% EMG Power
    fpass = [300,3000];
    trimmedEMG = RawData.data.EMG(1:min(analogExpectedLength,length(RawData.data.EMG)));
    [z,p,k] = butter(3,fpass/(ProcData.notes.analogSamplingRate/2));
    [sos,g] = zp2sos(z,p,k);
    filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
    kernelWidth = 0.5;
    smoothingKernel = gausswin(kernelWidth*ProcData.notes.analogSamplingRate)/sum(gausswin(kernelWidth*ProcData.notes.analogSamplingRate));
    EMGPwr = log10(conv(filtEMG.^2,smoothingKernel,'same'));
    resampEMG = resample(EMGPwr,ProcData.notes.dsFs,ProcData.notes.analogSamplingRate);
    ProcData.data.EMG.emg = resampEMG;
    %% EMG Signal
    fpass = [10,100];
    [z,p,k] = butter(3,fpass/(ProcData.notes.analogSamplingRate/2));
    [sos,g] = zp2sos(z,p,k);
    filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
    resampEMGSignal = resample(filtEMG,ProcData.notes.dsFs,ProcData.notes.analogSamplingRate);
    ProcData.data.EMG.emgSignal = resampEMGSignal;
    %% Pupil
    if isfield(RawData.data,'Pupil')==1
    % ProcData.data.Pupil.pupilArea = RawData.data.pupilArea; 
    % ProcData.data.Pupil.mmarea = RawData.data.pupilmmArea;
    ProcData.data.Pupil.Diameter = RawData.data.pupilDiameter;

        if isfield(RawData.data,'pupilMovement')==1
             ProcData.data.Pupil.Movement = RawData.data.pupilMovement;    
        end
    % ProcData.data.Pupil.mmDiameter = RawData.data.pupilmmDiameter;

    % ProcData.data.Pupil.Major = RawData.data.pupilpatchMajor;
    % ProcData.data.Pupil.Minor = RawData.data.pupilpatchMinor;
    % ProcData.data.Pupil.CentroidX = RawData.data.pupilpatchCentroidX;
    % ProcData.data.Pupil.CentroidY = RawData.data.pupilpatchCentroidY;

%     ProcData.data.Pupil.blinkFrames = RawData.data.Pupil.blinkFrames;
    % ProcData.data.Pupil.blinkInds = RawData.data.Pupil.blinkInds;
    % ProcData.data.Pupil.eyeROI = RawData.data.Pupil.eyeROI;
    % ProcData.data.Pupil.firstFrame = RawData.data.Pupil.firstFrame;
    end
    %% save the processed data
    save(procDataFile,'ProcData','-v7.3')
end

end
