function [] = CreateModelDataSet_FP_EEG(procDataFileIDs,NBins)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Arrange data into a table of most-relevant parameters for model training/classification
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    modelDataSetID = [procDataFileID(1:end-12) 'ModelData.mat'];
    load(procDataFileID)
    disp(['Creating model data set for ' procDataFileID '... (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')' ]); disp(' ')
    %% Create table to send into model
    variableNames = {'maxEEGDelta','maxEEGTheta','maxEEGAlpha','maxEEGBeta','maxEEGGamma',...
            'WhiskEvents','sumWhiskMotion',...
            'meanForce','stdForce',...
            'medEMG','stdEMG',...
            'meanMDiameter','minMDiameter',...
            'meanZDiameter','minZDiameter',...
            'sumEyeMotion','maxEyeMotion',...
            'meanCentroidX','maxCentroidX',...
            'meanCentroidY','maxCentroidY',...
            };%,...
  
    % pre-allocation
    maxEEGDelta_column = zeros(NBins,1);
    maxEEGBeta_column = zeros(NBins,1);
    maxEEGGamma_column = zeros(NBins,1);
    maxEEGTheta_column = zeros(NBins,1);
    maxEEGAlpha_column = zeros(NBins,1);

    numWhiskEvents_column = zeros(NBins,1);

    meanForce_column = zeros(NBins,1);
    stdForce_column = zeros(NBins,1);

    medEMG_column = zeros(NBins,1);
    stdEMG_column = zeros(NBins,1);

    avgMDiameterColumn = zeros(NBins,1);
    varMDiameterColumn = zeros(NBins,1);
    minMDiameterColumn = zeros(NBins,1);
    avgZDiameterColumn = zeros(NBins,1);
    varZDiameterColumn = zeros(NBins,1);
    minZDiameterColumn = zeros(NBins,1);
    sumEyeMotionColumn = zeros(NBins,1);
    maxEyeMotionColumn = zeros(NBins,1);
    varEyeMotionColumn = zeros(NBins,1);
    avgCentroidXColumn = zeros(NBins,1);
    varCentroidXColumn = zeros(NBins,1);
    maxCentroidXColumn = zeros(NBins,1);
    avgCentroidYColumn = zeros(NBins,1);
    varCentroidYColumn = zeros(NBins,1);
    maxCentroidYColumn = zeros(NBins,1);
    sumBinWhiskColumn = zeros(NBins,1);
    sumWhiskAngleColumn = zeros(NBins,1);
    varWhiskAngleColumn = zeros(NBins,1);

    % extract relevant parameters from each epoch
    for b = 1:length(maxEEGDelta_column)
        % EEG delta
        maxLHcortDelta = mean(cell2mat(ProcData.sleep.parameters.EEG_LH.specDeltaBandPower{b,1}));
        maxRHcortDelta = mean(cell2mat(ProcData.sleep.parameters.EEG_RH.specDeltaBandPower{b,1}));
        if maxLHcortDelta >= maxRHcortDelta
            maxEEGDelta_column(b,1) = maxLHcortDelta;
        else
            maxEEGDelta_column(b,1) = maxRHcortDelta;
        end
       
        % EEG theta
        maxLHcortTheta = mean(cell2mat(ProcData.sleep.parameters.EEG_LH.specThetaBandPower{b,1}));
        maxRHcortTheta = mean(cell2mat(ProcData.sleep.parameters.EEG_RH.specThetaBandPower{b,1}));
        if maxLHcortTheta >= maxRHcortTheta
            maxEEGTheta_column(b,1) = maxLHcortTheta;
        else
            maxEEGTheta_column(b,1) = maxRHcortTheta;
        end

        % EEG Alpha
        maxLHcortAlpha = mean(cell2mat(ProcData.sleep.parameters.EEG_LH.specAlphaBandPower{b,1}));
        maxRHcortAlpha = mean(cell2mat(ProcData.sleep.parameters.EEG_RH.specAlphaBandPower{b,1}));
        if maxLHcortAlpha >= maxRHcortAlpha
            maxEEGAlpha_column(b,1) = maxLHcortAlpha;
        else
            maxEEGAlpha_column(b,1) = maxRHcortAlpha;
        end

        % EEG beta
        maxLHcortBeta = mean(cell2mat(ProcData.sleep.parameters.EEG_LH.specBetaBandPower{b,1}));
        maxRHcortBeta = mean(cell2mat(ProcData.sleep.parameters.EEG_RH.specBetaBandPower{b,1}));
        if maxLHcortBeta >= maxRHcortBeta
            maxEEGBeta_column(b,1) = maxLHcortBeta;
        else
            maxEEGBeta_column(b,1) = maxRHcortBeta;
        end

        % EEG gamma
        maxLHcortGamma = mean(cell2mat(ProcData.sleep.parameters.EEG_LH.specGammaBandPower{b,1}));
        maxRHcortGamma = mean(cell2mat(ProcData.sleep.parameters.EEG_RH.specGammaBandPower{b,1}));
        if maxLHcortGamma >= maxRHcortGamma
            maxEEGGamma_column(b,1) = maxLHcortGamma;
        else
            maxEEGGamma_column(b,1) = maxRHcortGamma;
        end
      
        % number of binarized whisking events
        numWhiskEvents_column(b,1) = sum(ProcData.sleep.parameters.binWhiskerAngle{b,1});
        % number of binarized force sensor events
        meanForce_column(b,1) = mean(ProcData.sleep.parameters.ForceSensor{b,1});
        stdForce_column(b,1) = std(ProcData.sleep.parameters.ForceSensor{b,1});
        % average of the log of the EMG profile
        EMG = ProcData.sleep.parameters.EMG.power{b,1};
        medEMG_column(b,1) = median(EMG);
        stdEMG_column(b,1) = std(EMG);

        % pupil parameters
        avgMDiameterColumn(b,1) = mean(ProcData.sleep.parameters.Pupil.mmDiameter{b,1},'omitnan');
        varMDiameterColumn(b,1) = var(ProcData.sleep.parameters.Pupil.mmDiameter{b,1},'omitnan');
        minMDiameterColumn(b,1) = min(ProcData.sleep.parameters.Pupil.mmDiameter{b,1},[],'omitnan');
        avgZDiameterColumn(b,1) = mean(ProcData.sleep.parameters.Pupil.zDiameter{b,1},'omitnan');
        varZDiameterColumn(b,1) = var(ProcData.sleep.parameters.Pupil.zDiameter{b,1},'omitnan');
        minZDiameterColumn(b,1) = min(ProcData.sleep.parameters.Pupil.zDiameter{b,1},[],'omitnan');
        sumEyeMotionColumn(b,1) = sum(ProcData.sleep.parameters.Pupil.eyeMotion{b,1},'omitnan');
        maxEyeMotionColumn(b,1) = max(ProcData.sleep.parameters.Pupil.eyeMotion{b,1},[],'omitnan');
        varEyeMotionColumn(b,1) = var(ProcData.sleep.parameters.Pupil.eyeMotion{b,1},'omitnan');
        avgCentroidXColumn(b,1) = mean(ProcData.sleep.parameters.Pupil.CentroidX{b,1},'omitnan');
        varCentroidXColumn(b,1) = var(ProcData.sleep.parameters.Pupil.CentroidX{b,1},'omitnan');
        maxCentroidXColumn(b,1) = max(ProcData.sleep.parameters.Pupil.CentroidX{b,1},[],'omitnan');
        avgCentroidYColumn(b,1) = mean(ProcData.sleep.parameters.Pupil.CentroidY{b,1},'omitnan');
        varCentroidYColumn(b,1) = var(ProcData.sleep.parameters.Pupil.CentroidY{b,1},'omitnan');
        maxCentroidYColumn(b,1) = max(ProcData.sleep.parameters.Pupil.CentroidY{b,1},[],'omitnan');
        sumBinWhiskColumn(b,1) = sum(ProcData.sleep.parameters.binWhiskerAngle{b,1},'omitnan');
        sumWhiskAngleColumn(b,1) = sum(ProcData.sleep.parameters.Pupil.whiskerMotion{b,1},'omitnan');
        varWhiskAngleColumn(b,1) = var(ProcData.sleep.parameters.Pupil.whiskerMotion{b,1},'omitnan');
    end
    % create table
    paramsTable = table(maxEEGDelta_column,maxEEGTheta_column,maxEEGAlpha_column,maxEEGBeta_column,maxEEGGamma_column,...
            sumBinWhiskColumn,sumWhiskAngleColumn,...
            meanForce_column,stdForce_column,...
            medEMG_column,stdEMG_column,...
            avgMDiameterColumn,minMDiameterColumn,...
            avgZDiameterColumn,minZDiameterColumn,...
            sumEyeMotionColumn,maxEyeMotionColumn,...
            avgCentroidXColumn,maxCentroidXColumn,...
            avgCentroidYColumn,maxCentroidYColumn,...
            'VariableNames',variableNames);
    save(modelDataSetID,'paramsTable')
end

end
