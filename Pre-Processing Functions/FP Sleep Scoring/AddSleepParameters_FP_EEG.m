function [] = AddSleepParameters_FP_EEG(procDataFileIDs,RestingBaselines,baselineType,NBins)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Organize data into appropriate bins for sleep scoring characterization
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Adding sleep scoring parameters to ' procDataFileID '... (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
    [~, fileDate, ~] = GetFileInfo_FP(procDataFileID);
    strDay = ConvertDate_FP(fileDate);
    load(procDataFileID)
%     ProcData.sleep = rmfield(ProcData.sleep,'parameters');
    specDataFileID = [procDataFileID(1:end-12) 'SpecDataB.mat'];
    load(specDataFileID)
    %% BLOCK PURPOSE: Create folder for the Neural data of each electrode

    % EEG delta
    LH_Delta = ProcData.data.EEG_LH.deltaBandPower;
    RH_Delta = ProcData.data.EEG_RH.deltaBandPower;
    LH_baselineDelta = RestingBaselines.(baselineType).EEG_LH.deltaBandPower.(strDay).mean;
    RH_baselineDelta = RestingBaselines.(baselineType).EEG_RH.deltaBandPower.(strDay).mean;
    LH_DeltaNeuro = (LH_Delta-LH_baselineDelta)/LH_baselineDelta;
    RH_DeltaNeuro = (RH_Delta-RH_baselineDelta)/RH_baselineDelta;  

    % EEG theta
    LH_Theta = ProcData.data.EEG_LH.thetaBandPower;
    RH_Theta = ProcData.data.EEG_RH.thetaBandPower;
    LH_baselineTheta = RestingBaselines.(baselineType).EEG_LH.thetaBandPower.(strDay).mean;
    RH_baselineTheta = RestingBaselines.(baselineType).EEG_LH.thetaBandPower.(strDay).mean;
    LH_ThetaNeuro = (LH_Theta-LH_baselineTheta)/LH_baselineTheta;
    RH_ThetaNeuro = (RH_Theta-RH_baselineTheta)/RH_baselineTheta;

    % EEG alpha
    LH_Alpha = ProcData.data.EEG_LH.alphaBandPower;
    RH_Alpha = ProcData.data.EEG_RH.alphaBandPower;
    LH_baselineAlpha = RestingBaselines.(baselineType).EEG_LH.alphaBandPower.(strDay).mean;
    RH_baselineAlpha = RestingBaselines.(baselineType).EEG_LH.alphaBandPower.(strDay).mean;
    LH_AlphaNeuro = (LH_Alpha-LH_baselineAlpha)/LH_baselineAlpha;
    RH_AlphaNeuro = (RH_Alpha-RH_baselineAlpha)/RH_baselineAlpha;   
    % EEG beta
    LH_Beta = ProcData.data.EEG_LH.betaBandPower;
    RH_Beta = ProcData.data.EEG_RH.betaBandPower;
    LH_baselineBeta = RestingBaselines.(baselineType).EEG_LH.betaBandPower.(strDay).mean;
    RH_baselineBeta = RestingBaselines.(baselineType).EEG_LH.betaBandPower.(strDay).mean;
    LH_BetaNeuro = (LH_Beta-LH_baselineBeta)/LH_baselineBeta;
    RH_BetaNeuro = (RH_Beta-RH_baselineBeta)/RH_baselineBeta;

    % EEG gamma
    LH_Gamma = ProcData.data.EEG_LH.gammaBandPower;
    RH_Gamma = ProcData.data.EEG_RH.gammaBandPower;
    LH_baselineGamma = RestingBaselines.(baselineType).EEG_LH.gammaBandPower.(strDay).mean;
    RH_baselineGamma = RestingBaselines.(baselineType).EEG_LH.gammaBandPower.(strDay).mean;
    LH_GammaNeuro = (LH_Gamma-LH_baselineGamma)/LH_baselineGamma;
    RH_GammaNeuro = (RH_Gamma-RH_baselineGamma)/RH_baselineGamma;

    % EEG EEG
    LH_EEG = ProcData.data.EEG_LH.EEGPower;
    RH_EEG = ProcData.data.EEG_RH.EEGPower;
    LH_baselineEEG = RestingBaselines.(baselineType).EEG_LH.EEGPower.(strDay).mean;
    RH_baselineEEG = RestingBaselines.(baselineType).EEG_LH.EEGPower.(strDay).mean;
    LH_EEGNeuro = (LH_EEG-LH_baselineEEG)/LH_baselineEEG;
    RH_EEGNeuro = (RH_EEG-RH_baselineEEG)/RH_baselineEEG;
    % Divide the neural signals into five second bins and put them in a cell array

    LH_tempDeltaStruct = cell(NBins,1);
    RH_tempDeltaStruct = cell(NBins,1);
    LH_tempThetaStruct = cell(NBins,1);
    RH_tempThetaStruct = cell(NBins,1);
    LH_tempAlphaStruct = cell(NBins,1);
    RH_tempAlphaStruct = cell(NBins,1);
    LH_tempBetaStruct = cell(NBins,1);
    RH_tempBetaStruct = cell(NBins,1);
    LH_tempGammaStruct = cell(NBins,1);
    RH_tempGammaStruct = cell(NBins,1);
    LH_tempEEGStruct = cell(NBins,1);
    RH_tempEEGStruct = cell(NBins,1);
    % loop through all samples across the full duration in 5 second bins (NBins total)
    for b = 1:NBins
        if b == 1
            % EEG
            LH_tempDeltaStruct(b,1) = {LH_DeltaNeuro(b:150)};
            RH_tempDeltaStruct(b,1) = {RH_DeltaNeuro(b:150)};
            LH_tempThetaStruct(b,1) = {LH_ThetaNeuro(b:150)};
            RH_tempThetaStruct(b,1) = {RH_ThetaNeuro(b:150)};
            LH_tempAlphaStruct(b,1) = {LH_AlphaNeuro(b:150)};
            RH_tempAlphaStruct(b,1) = {RH_AlphaNeuro(b:150)};
            LH_tempBetaStruct(b,1) = {LH_BetaNeuro(b:150)};
            RH_tempBetaStruct(b,1) = {RH_BetaNeuro(b:150)};
            LH_tempGammaStruct(b,1) = {LH_GammaNeuro(b:150)};
            RH_tempGammaStruct(b,1) = {RH_GammaNeuro(b:150)};
            LH_tempEEGStruct(b,1) = {LH_EEGNeuro(b:150)};
            RH_tempEEGStruct(b,1) = {RH_EEGNeuro(b:150)};
        elseif b == NBins
            % EEG
            LH_tempDeltaStruct(b,1) = {LH_DeltaNeuro((((150*(b - 1)) + 1)):end)};
            RH_tempDeltaStruct(b,1) = {RH_DeltaNeuro((((150*(b - 1)) + 1)):end)};
            LH_tempThetaStruct(b,1) = {LH_ThetaNeuro((((150*(b - 1)) + 1)):end)};
            RH_tempThetaStruct(b,1) = {RH_ThetaNeuro((((150*(b - 1)) + 1)):end)};
            LH_tempAlphaStruct(b,1) = {LH_AlphaNeuro((((150*(b - 1)) + 1)):end)};
            RH_tempAlphaStruct(b,1) = {RH_AlphaNeuro((((150*(b - 1)) + 1)):end)};
            LH_tempBetaStruct(b,1) = {LH_BetaNeuro((((150*(b - 1)) + 1)):end)};
            RH_tempBetaStruct(b,1) = {RH_BetaNeuro((((150*(b - 1)) + 1)):end)};
            LH_tempGammaStruct(b,1) = {LH_GammaNeuro((((150*(b - 1)) + 1)):end)};
            RH_tempGammaStruct(b,1) = {RH_GammaNeuro((((150*(b - 1)) + 1)):end)};
            LH_tempEEGStruct(b,1) = {LH_EEGNeuro((((150*(b - 1)) + 1)):end)};
            RH_tempEEGStruct(b,1) = {RH_EEGNeuro((((150*(b - 1)) + 1)):end)};
        else
            % EEG
            LH_tempDeltaStruct(b,1) = {LH_DeltaNeuro((((150*(b - 1)) + 1)):(150*b))};
            RH_tempDeltaStruct(b,1) = {RH_DeltaNeuro((((150*(b - 1)) + 1)):(150*b))};
            LH_tempThetaStruct(b,1) = {LH_ThetaNeuro((((150*(b - 1)) + 1)):(150*b))};
            RH_tempThetaStruct(b,1) = {RH_ThetaNeuro((((150*(b - 1)) + 1)):(150*b))};
            LH_tempAlphaStruct(b,1) = {LH_AlphaNeuro((((150*(b - 1)) + 1)):(150*b))};
            RH_tempAlphaStruct(b,1) = {RH_AlphaNeuro((((150*(b - 1)) + 1)):(150*b))};
            LH_tempBetaStruct(b,1) = {LH_BetaNeuro((((150*(b - 1)) + 1)):(150*b))};
            RH_tempBetaStruct(b,1) = {RH_BetaNeuro((((150*(b - 1)) + 1)):(150*b))};
            LH_tempGammaStruct(b,1) = {LH_GammaNeuro((((150*(b - 1)) + 1)):(150*b))};
            RH_tempGammaStruct(b,1) = {RH_GammaNeuro((((150*(b - 1)) + 1)):(150*b))};
            LH_tempEEGStruct(b,1) = {LH_EEGNeuro((((150*(b - 1)) + 1)):(150*b))};
            RH_tempEEGStruct(b,1) = {RH_EEGNeuro((((150*(b - 1)) + 1)):(150*b))};
        end
    end
  
    % save EEG data under ProcData file
    ProcData.sleep.parameters.EEG_LH.deltaBandPower = LH_tempDeltaStruct;
    ProcData.sleep.parameters.EEG_RH.deltaBandPower = RH_tempDeltaStruct;
    ProcData.sleep.parameters.EEG_LH.thetaBandPower = LH_tempThetaStruct;
    ProcData.sleep.parameters.EEG_RH.thetaBandPower = RH_tempThetaStruct;
    ProcData.sleep.parameters.EEG_LH.alphaBandPower = LH_tempAlphaStruct;
    ProcData.sleep.parameters.EEG_RH.alphaBandPower = RH_tempAlphaStruct;
    ProcData.sleep.parameters.EEG_LH.betaBandPower = LH_tempBetaStruct;
    ProcData.sleep.parameters.EEG_RH.betaBandPower = RH_tempBetaStruct;
    ProcData.sleep.parameters.EEG_LH.gammaBandPower = LH_tempGammaStruct;
    ProcData.sleep.parameters.EEG_RH.gammaBandPower = RH_tempGammaStruct;
    ProcData.sleep.parameters.EEG_LH.EEGPower = LH_tempEEGStruct;
    ProcData.sleep.parameters.EEG_RH.EEGPower = RH_tempEEGStruct;
    %% BLOCK PURPOSE: Create folder for the Neural spectrogram data of each electrode
    trialDuration_sec = ProcData.notes.trialDuration_sec;   % sec
    offset = 2.5;   % sec
    binWidth = 5;   % sec
    T = round(SpecData.EEG_LH.T,1);
    F = SpecData.EEG_LH.F;
    specLH = SpecData.EEG_LH.normS;
    specRH = SpecData.EEG_RH.normS;
    freqFloor = floor(F);
    % delta
    deltaLow = freqFloor == 1;%1;
    deltaHigh = freqFloor == 4;
    deltaLowStart = find(deltaLow,1,'first');
    deltaLowEnd = find(deltaHigh,1,'last');
    deltaSpecLH = specLH(deltaLowStart:deltaLowEnd,:);
    deltaSpecRH = specRH(deltaLowStart:deltaLowEnd,:);
    meanDeltaSpecLH = mean(deltaSpecLH,1);
    meanDeltaSpecRH = mean(deltaSpecRH,1);
    % theta
    thetaLow = freqFloor == 4;
    thetaHigh = freqFloor == 9;%10;
    thetaLowStart = find(thetaLow,1,'first');
    thetaLowEnd = find(thetaHigh,1,'last');
    thetaSpecLH = specLH(thetaLowStart:thetaLowEnd,:);
    thetaSpecRH = specRH(thetaLowStart:thetaLowEnd,:);
    meanThetaSpecLH = mean(thetaSpecLH,1);
    meanThetaSpecRH = mean(thetaSpecRH,1);
    % alpha
    alphaLow = freqFloor == 10;
    alphaHigh = freqFloor == 13;
    alphaLowStart = find(alphaLow,1,'first');
    alphaLowEnd = find(alphaHigh,1,'last');
    alphaSpecLH = specLH(alphaLowStart:alphaLowEnd,:);
    alphaSpecRH = specRH(alphaLowStart:alphaLowEnd,:);
    meanAlphaSpecLH = mean(alphaSpecLH,1);
    meanAlphaSpecRH = mean(alphaSpecRH,1);
    % beta
    betaLow = freqFloor == 13;
    betaHigh = freqFloor == 30;
    betaLowStart = find(betaLow,1,'first');
    betaLowEnd = find(betaHigh,1,'last');
    betaSpecLH = specLH(betaLowStart:betaLowEnd,:);
    betaSpecRH = specRH(betaLowStart:betaLowEnd,:);
    meanBetaSpecLH = mean(betaSpecLH,1);
    meanBetaSpecRH = mean(betaSpecRH,1);
    % gamma
    gammaLow = freqFloor == 30;
    gammaHigh = freqFloor == 99;
    gammaLowStart = find(gammaLow,1,'first');
    gammaLowEnd = find(gammaHigh,1,'last');
    gammaSpecLH = specLH(gammaLowStart:gammaLowEnd,:);
    gammaSpecRH = specRH(gammaLowStart:gammaLowEnd,:);
    meanGammaSpecRH = mean(gammaSpecRH,1);
    meanGammaSpecLH = mean(gammaSpecLH,1);

    LH_tempDeltaSpecStruct = cell(NBins,1);
    RH_tempDeltaSpecStruct = cell(NBins,1);
    LH_tempThetaSpecStruct = cell(NBins,1);
    RH_tempThetaSpecStruct = cell(NBins,1);
    LH_tempAlphaSpecStruct = cell(NBins,1);
    RH_tempAlphaSpecStruct = cell(NBins,1);
    LH_tempBetaSpecStruct = cell(NBins,1);
    RH_tempBetaSpecStruct = cell(NBins,1);
    LH_tempGammaSpecStruct = cell(NBins,1);
    RH_tempGammaSpecStruct = cell(NBins,1);
    % loop through all samples across the 15 minutes in 5 second bins (NBins total)
    for c = 1:NBins
        if c == 1
            startTime = offset;
            startTime_index = find(T == startTime);
            endTime = 5;
            [~,endTime_index] = min(abs(T - endTime));
            % EEG
            LH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecLH(startTime_index:endTime_index)};
            RH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecRH(startTime_index:endTime_index)};
            LH_tempThetaSpecStruct{c,1} = {meanThetaSpecLH(startTime_index:endTime_index)};
            RH_tempThetaSpecStruct{c,1} = {meanThetaSpecRH(startTime_index:endTime_index)};
            LH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecLH(startTime_index:endTime_index)};
            RH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecRH(startTime_index:endTime_index)};
            LH_tempBetaSpecStruct{c,1} = {meanBetaSpecLH(startTime_index:endTime_index)};
            RH_tempBetaSpecStruct{c,1} = {meanBetaSpecRH(startTime_index:endTime_index)};
            LH_tempGammaSpecStruct{c,1} = {meanGammaSpecLH(startTime_index:endTime_index)};
            RH_tempGammaSpecStruct{c,1} = {meanGammaSpecRH(startTime_index:endTime_index)};
        elseif c == NBins
            startTime = trialDuration_sec - 5;
            [~,startTime_index] = min(abs(T - startTime));
            endTime = trialDuration_sec - offset;
            [~,endTime_index] = min(abs(T - endTime));
            % EEG
            LH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecLH(startTime_index:endTime_index)};
            RH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecRH(startTime_index:endTime_index)};
            LH_tempThetaSpecStruct{c,1} = {meanThetaSpecLH(startTime_index:endTime_index)};
            RH_tempThetaSpecStruct{c,1} = {meanThetaSpecRH(startTime_index:endTime_index)};
            LH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecLH(startTime_index:endTime_index)};
            RH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecRH(startTime_index:endTime_index)};
            LH_tempBetaSpecStruct{c,1} = {meanBetaSpecLH(startTime_index:endTime_index)};
            RH_tempBetaSpecStruct{c,1} = {meanBetaSpecRH(startTime_index:endTime_index)};
            LH_tempGammaSpecStruct{c,1} = {meanGammaSpecLH(startTime_index:endTime_index)};
            RH_tempGammaSpecStruct{c,1} = {meanGammaSpecRH(startTime_index:endTime_index)};
        else
            startTime = binWidth*(c - 1);
            [~,startTime_index] = min(abs(T - startTime));
            endTime = binWidth*c;
            [~,endTime_index] = min(abs(T - endTime));
            % EEG
            LH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecLH(startTime_index + 1:endTime_index + 1)};
            RH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecRH(startTime_index + 1:endTime_index + 1)};
            LH_tempThetaSpecStruct{c,1} = {meanThetaSpecLH(startTime_index + 1:endTime_index + 1)};
            RH_tempThetaSpecStruct{c,1} = {meanThetaSpecRH(startTime_index + 1:endTime_index + 1)};
            LH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecLH(startTime_index + 1:endTime_index + 1)};
            RH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecRH(startTime_index + 1:endTime_index + 1)};
            LH_tempBetaSpecStruct{c,1} = {meanBetaSpecLH(startTime_index + 1:endTime_index + 1)};
            RH_tempBetaSpecStruct{c,1} = {meanBetaSpecRH(startTime_index + 1:endTime_index + 1)};
            LH_tempGammaSpecStruct{c,1} = {meanGammaSpecLH(startTime_index + 1:endTime_index + 1)};
            RH_tempGammaSpecStruct{c,1} = {meanGammaSpecRH(startTime_index + 1:endTime_index + 1)};
        end
    end

    % save EEG data under ProcData file
    ProcData.sleep.parameters.EEG_LH.specDeltaBandPower = LH_tempDeltaSpecStruct;
    ProcData.sleep.parameters.EEG_RH.specDeltaBandPower = RH_tempDeltaSpecStruct;
    ProcData.sleep.parameters.EEG_LH.specThetaBandPower = LH_tempThetaSpecStruct;
    ProcData.sleep.parameters.EEG_RH.specThetaBandPower = RH_tempThetaSpecStruct;
    ProcData.sleep.parameters.EEG_LH.specAlphaBandPower = LH_tempAlphaSpecStruct;
    ProcData.sleep.parameters.EEG_RH.specAlphaBandPower = RH_tempAlphaSpecStruct;
    ProcData.sleep.parameters.EEG_LH.specBetaBandPower = LH_tempBetaSpecStruct;
    ProcData.sleep.parameters.EEG_RH.specBetaBandPower = RH_tempBetaSpecStruct;
    ProcData.sleep.parameters.EEG_LH.specGammaBandPower = LH_tempGammaSpecStruct;
    ProcData.sleep.parameters.EEG_RH.specGammaBandPower = RH_tempGammaSpecStruct;
    
    %% BLOCK PURPOSE: Create folder for binarized whisking and binarized force sensor
    binWhiskerAngle = ProcData.data.binWhiskerAngle;
    binForceSensor = ProcData.data.binForceSensor;
    ForceSensor = ProcData.data.forceSensor;
    whiskerAngle = ProcData.data.whiskerAngle;
    whiskerAcceleration = diff(whiskerAngle,2);
    % Find the number of whiskerBins due to frame drops.
    whiskerBinNumber = NBins;
    % Divide the signal into five second bins and put them in a cell array
    tempWhiskerStruct = cell(whiskerBinNumber,1);
    tempWhiskerAccelStruct = cell(whiskerBinNumber,1);
    tempBinWhiskerStruct = cell(whiskerBinNumber,1);
    tempForceStruct = cell(whiskerBinNumber,1);
    tempForceRStruct = cell(whiskerBinNumber,1);
    for whiskerBins = 1:whiskerBinNumber
        if whiskerBins == 1
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle(whiskerBins:150)};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration(whiskerBins:150)};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle(whiskerBins:150)};
            tempForceStruct(whiskerBins, 1) = {binForceSensor(whiskerBins:150)};
            tempForceRStruct(whiskerBins, 1) = {ForceSensor(whiskerBins:150)};
        elseif whiskerBins == whiskerBinNumber
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle((((150*(whiskerBins-1)) + 1)):end)};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration((((150*(whiskerBins-1)) + 1)):end)};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle((((150*(whiskerBins-1)) + 1)):end)};
            tempForceStruct(whiskerBins, 1) = {binForceSensor((((150*(whiskerBins-1)) + 1)):end)};
            tempForceRStruct(whiskerBins, 1) = {ForceSensor((((150*(whiskerBins-1)) + 1)):end)};
        else
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
            tempForceStruct(whiskerBins, 1) = {binForceSensor((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
            tempForceRStruct(whiskerBins, 1) = {ForceSensor((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
        end
    end
    % save whisker and force sensor data under ProcData file
    ProcData.sleep.parameters.whiskerAngle = tempWhiskerStruct;
    ProcData.sleep.parameters.whiskerAcceleration = tempWhiskerAccelStruct;
    ProcData.sleep.parameters.binWhiskerAngle = tempBinWhiskerStruct;
    ProcData.sleep.parameters.binForceSensor = tempForceStruct;
    ProcData.sleep.parameters.ForceSensor = tempForceRStruct;
        %% add pupil parameters
    ProcData.data.Pupil.mmPerPixel = 0.018; % used to convert pixels to mm
        if strcmp(ProcData.data.Pupil.newDiameterCheck,'y') == true
        % create fields for the data
        dataTypes = {'pupilArea','Diameter','mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
            for aa = 1:length(dataTypes)
                dataType = dataTypes{1,aa};
                samplingRate = ProcData.notes.dsFs;
                [z,p,k] = butter(4,1/(samplingRate/2),'low');
                [sos,g] = zp2sos(z,p,k);
                if strcmp(dataType,'eyeMotion') == true
                    [z,p,k] = butter(4,10/(samplingRate/2),'low');
                    [sos,g] = zp2sos(z,p,k);
                    try
                        centroidX = filtfilt(sos,g,ProcData.data.Pupil.CentroidX - RestingBaselines.manualSelection.Pupil.CentroidX.(strDay).mean);
                        centroidY = filtfilt(sos,g,ProcData.data.Pupil.CentroidY - RestingBaselines.manualSelection.Pupil.CentroidY.(strDay).mean);
                    catch
                        centroidX = ProcData.data.Pupil.CentroidX - RestingBaselines.manualSelection.Pupil.CentroidX.(strDay).mean;
                        centroidY = ProcData.data.Pupil.CentroidY - RestingBaselines.manualSelection.Pupil.CentroidY.(strDay).mean;
                    end
                    for dd = 1:length(centroidX)
                        if dd == 1
                            distanceTraveled(1,dd) = 0;
                        else
                            xPt = [centroidX(1,dd - 1),centroidY(1,dd - 1)];
                            yPt = [centroidX(1,dd),centroidY(1,dd)];
                            distanceTraveled(1,dd) = pdist([xPt;yPt],'euclidean');
                        end
                    end

                    data.(dataType).data = distanceTraveled*ProcData.data.Pupil.mmPerPixel;
                    ProcData.data.Pupil.distanceTraveled = distanceTraveled*ProcData.data.Pupil.mmPerPixel;
                elseif strcmp(dataType,'CentroidX') == true
                    try
                        data.(dataType).data = filtfilt(sos,g,ProcData.data.Pupil.CentroidX - RestingBaselines.manualSelection.Pupil.CentroidX.(strDay).mean);
                    catch
                        data.(dataType).data = ProcData.data.Pupil.CentroidX - RestingBaselines.manualSelection.Pupil.CentroidX.(strDay).mean;
                    end
                elseif strcmp(dataType,'CentroidY') == true
                    try
                        data.(dataType).data = filtfilt(sos,g,ProcData.data.Pupil.CentroidY - RestingBaselines.manualSelection.Pupil.CentroidY.(strDay).mean);
                    catch
                        data.(dataType).data = ProcData.data.Pupil.CentroidY - RestingBaselines.manualSelection.Pupil.CentroidY.(strDay).mean;
                    end
                elseif strcmp(dataType,'whiskerMotion')
                    data.(dataType).data = ProcData.data.whiskerAngle.^2;
                else
                    try
                        data.(dataType).data = filtfilt(sos,g,ProcData.data.Pupil.(dataType));
                    catch
                        data.(dataType).data = ProcData.data.Pupil.(dataType);
                    end
                end
                data.(dataType).struct = cell(NBins,1);
                % loop through all samples across the 15 minutes in 5 second bins (NBins total)
                for b = 1:NBins
                    if b == 1
                        data.(dataType).struct(b,1) = {data.(dataType).data(b:150)};
                    elseif b == NBins
                        data.(dataType).struct(b,1) = {data.(dataType).data((((150*(b - 1)) + 1)):(150*b))};%end)};
                    else
                        data.(dataType).struct(b,1) = {data.(dataType).data((((150*(b - 1)) + 1)):(150*b))};
                    end
                end
                ProcData.sleep.parameters.Pupil.(dataType) = data.(dataType).struct;
            end
        end
    %% Create folder for the EMG
    EMG = ProcData.data.EMG.emg;
    normEMG = EMG - RestingBaselines.(baselineType).EMG.power.(strDay).mean;  
    tempEMGStruct = cell(NBins,1);
    for EMGBins = 1:NBins
        if EMGBins == 1
            tempEMGStruct(EMGBins,1) = {normEMG(EMGBins:150)};
        else
            tempEMGStruct(EMGBins,1) = {normEMG((((150*(EMGBins-1)) + 1)):(150*EMGBins))};
        end
    end
    % save EMG data under ProcData file
    ProcData.sleep.parameters.EMG.emg = tempEMGStruct;

    EMG = ProcData.data.EMG.power;
    normEMG = EMG - RestingBaselines.(baselineType).EMG.emg.(strDay).mean;  
    tempEMGStruct = cell(NBins,1);
    for EMGBins = 1:NBins
        if EMGBins == 1
            tempEMGStruct(EMGBins,1) = {normEMG(EMGBins:150)};
        else
            tempEMGStruct(EMGBins,1) = {normEMG((((150*(EMGBins-1)) + 1)):(150*EMGBins))};
        end
    end
    % save EMG data under ProcData file
    ProcData.sleep.parameters.EMG.power = tempEMGStruct;      
    %% BLOCK PURPOSE: Create folder for the left and right Rhodamine data
%     LH_GFP = ProcData.data.GFP.LH;
    RH_GFP = ProcData.data.GFP.RH; 
%     LH_Rhodamine = ProcData.data.Rhodamine.LH;
    RH_Rhodamine = ProcData.data.Rhodamine.RH;

%     LH_NormRhodamine = (LH_Rhodamine - RestingBaselines.(baselineType).Rhodamine.LH.(strDay).mean)/RestingBaselines.(baselineType).Rhodamine.LH.(strDay).std;
    RH_NormRhodamine = (RH_Rhodamine - RestingBaselines.(baselineType).Rhodamine.RH.(strDay).mean)/RestingBaselines.(baselineType).Rhodamine.RH.(strDay).std;
%     LH_NormGFP = (LH_GFP - RestingBaselines.(baselineType).GFP.LH.(strDay).mean)/RestingBaselines.(baselineType).GFP.LH.(strDay).std;
    RH_NormGFP = (RH_GFP - RestingBaselines.(baselineType).GFP.RH.(strDay).mean)/RestingBaselines.(baselineType).GFP.RH.(strDay).std;
    

%     LH_tempGFPStruct = cell(NBins,1);
    RH_tempGFPStruct = cell(NBins,1);  
%     RhodamineLH_tempRhodamineStruct = cell(NBins,1);
    RhodamineRH_tempRhodamineStruct = cell(NBins,1);

    for RhodamineBins = 1:NBins
        if RhodamineBins == 1
%             LH_tempGFPStruct(RhodamineBins,1) = {LH_NormGFP(RhodamineBins:150)};  % Samples 1 to 150
            RH_tempGFPStruct(RhodamineBins,1) = {RH_NormGFP(RhodamineBins:150)};

%             RhodamineLH_tempRhodamineStruct(RhodamineBins,1) = {LH_NormRhodamine(RhodamineBins:150)};  % Samples 1 to 150
            RhodamineRH_tempRhodamineStruct(RhodamineBins,1) = {RH_NormRhodamine(RhodamineBins:150)};
        else
%             LH_tempGFPStruct(RhodamineBins,1) = {LH_NormGFP((((150*(RhodamineBins-1)) + 1)):(150*RhodamineBins))};  % Samples 151 to 300, etc...
            RH_tempGFPStruct(RhodamineBins,1) = {RH_NormGFP((((150*(RhodamineBins-1)) + 1)):(150*RhodamineBins))};
%             RhodamineLH_tempRhodamineStruct(RhodamineBins,1) = {LH_NormRhodamine((((150*(RhodamineBins-1)) + 1)):(150*RhodamineBins))};  % Samples 151 to 300, etc...
            RhodamineRH_tempRhodamineStruct(RhodamineBins,1) = {RH_NormRhodamine((((150*(RhodamineBins-1)) + 1)):(150*RhodamineBins))};

        end
    end
    % save hemodynamic data under ProcData file
%     ProcData.sleep.parameters.GFP.LH = LH_tempGFPStruct;
    ProcData.sleep.parameters.GFP.RH = RH_tempGFPStruct;
%     ProcData.sleep.parameters.Rhodamine.LH = RhodamineLH_tempRhodamineStruct;
    ProcData.sleep.parameters.Rhodamine.RH = RhodamineRH_tempRhodamineStruct;
    save(procDataFileID,'ProcData');
end

end
