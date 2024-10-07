function [] = CreateModelDataSet_FP_Spec(procDataFileIDs,NBins)
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
    % pre-allocation
    CortDeltaLH_column = zeros(NBins,26);
    CortThetaLH_column = zeros(NBins,26);
%     CortAlphaLH_column = zeros(NBins,1);
    CortBetaLH_column = zeros(NBins,26);
    CortGammaLH_column = zeros(NBins,26);

    CortDeltaRH_column = zeros(NBins,26);
    CortThetaRH_column = zeros(NBins,26);
%     CortAlphaRH_column = zeros(NBins,1);
    CortBetaRH_column = zeros(NBins,26);
    CortGammaRH_column = zeros(NBins,26);

    HippDelta_column = zeros(NBins,26);
    HippTheta_column = zeros(NBins,26);
    HippGamma_column = zeros(NBins,26);
    WhiskEvents_column = zeros(NBins,26);
    ForceEvents_column = zeros(NBins,26);
    EMG_column = zeros(NBins,26);
%     avgHeartRate_column = zeros(NBins,1);
    % extract relevant parameters from each epoch
    for b = 1:length(CortDeltaLH_column)
        if b==1 || b == length(CortDeltaLH_column)
            %%
            % Cortical delta
            CortDeltaLH_column(b,:) = resample(cell2mat(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b}),26, length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
            CortDeltaRH_column(b,:) = resample(cell2mat(ProcData.sleep.parameters.cortical_RH.specDeltaBandPower{b}),26,length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
            % Cortical theta
            CortThetaLH_column(b,:) = resample(cell2mat(ProcData.sleep.parameters.cortical_LH.specThetaBandPower{b}),26,length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
            CortThetaRH_column(b,:) = resample(cell2mat(ProcData.sleep.parameters.cortical_RH.specThetaBandPower{b}),26,length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
            % Cortical alpha
    %         CortAlphaLH_column(b,1) = resample(cell2mat(ProcData.sleep.parameters.cortical_LH.specAlphaBandPower{b}),26,length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
    %         CortAlphaRH_column(b,1) = resample(cell2mat(ProcData.sleep.parameters.cortical_RH.specAlphaBandPower{b}),26,length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
            % Cortical beta
            CortBetaLH_column(b,:) = resample(cell2mat(ProcData.sleep.parameters.cortical_LH.specBetaBandPower{b}),26,length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
            CortBetaRH_column(b,:) = resample(cell2mat(ProcData.sleep.parameters.cortical_RH.specBetaBandPower{b}),26,length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
            % Cortical gamma
            CortGammaLH_column(b,:) = resample(cell2mat(ProcData.sleep.parameters.cortical_LH.specGammaBandPower{b}),26,length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
            CortGammaRH_column(b,:) = resample(cell2mat(ProcData.sleep.parameters.cortical_RH.specGammaBandPower{b}),26,length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
            % Hippocampal delta
            HippDelta_column(b,:) = resample(cell2mat(ProcData.sleep.parameters.hippocampus.specDeltaBandPower{b}),26,length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
            % Hippocampal theta
            HippTheta_column(b,:) = resample(cell2mat(ProcData.sleep.parameters.hippocampus.specThetaBandPower{b}),26,length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
            % Hippocampal theta
            HippGamma_column(b,:) = resample(cell2mat(ProcData.sleep.parameters.hippocampus.specGammaBandPower{b}),26,length(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b, 1}{1}));
        else
            % Cortical delta
            CortDeltaLH_column(b,:) = cell2mat(ProcData.sleep.parameters.cortical_LH.specDeltaBandPower{b});
            CortDeltaRH_column(b,:) = cell2mat(ProcData.sleep.parameters.cortical_RH.specDeltaBandPower{b});
            % Cortical theta
            CortThetaLH_column(b,:) = cell2mat(ProcData.sleep.parameters.cortical_LH.specThetaBandPower{b});
            CortThetaRH_column(b,:) = cell2mat(ProcData.sleep.parameters.cortical_RH.specThetaBandPower{b});
            % Cortical alpha
    %         CortAlphaLH_column(b,1) = cell2mat(ProcData.sleep.parameters.cortical_LH.specAlphaBandPower{b});
    %         CortAlphaRH_column(b,1) = cell2mat(ProcData.sleep.parameters.cortical_RH.specAlphaBandPower{b});
            % Cortical beta
            CortBetaLH_column(b,:) = cell2mat(ProcData.sleep.parameters.cortical_LH.specBetaBandPower{b});
            CortBetaRH_column(b,:) = cell2mat(ProcData.sleep.parameters.cortical_RH.specBetaBandPower{b});
            % Cortical gamma
            CortGammaLH_column(b,:) = cell2mat(ProcData.sleep.parameters.cortical_LH.specGammaBandPower{b});
            CortGammaRH_column(b,:) = cell2mat(ProcData.sleep.parameters.cortical_RH.specGammaBandPower{b});
            % Hippocampal delta
            HippDelta_column(b,:) = cell2mat(ProcData.sleep.parameters.hippocampus.specDeltaBandPower{b});
            % Hippocampal theta
            HippTheta_column(b,:) = cell2mat(ProcData.sleep.parameters.hippocampus.specThetaBandPower{b});
            % Hippocampal theta
            HippGamma_column(b,:) = cell2mat(ProcData.sleep.parameters.hippocampus.specGammaBandPower{b});
        end
        %%
        if b==length(CortDeltaLH_column)
            % number of binarized whisking events
            WhiskEvents= double([(ProcData.sleep.parameters.binWhiskerAngle{b}),0,0]);
            % number of binarized force sensor events
            ForceEvents = double([(ProcData.sleep.parameters.binForceSensor{b}),0]);  
        else
            % number of binarized whisking events
            WhiskEvents = double([(ProcData.sleep.parameters.binWhiskerAngle{b})]);
            % number of binarized force sensor events
            ForceEvents = double([(ProcData.sleep.parameters.binForceSensor{b})]);
        end
        % average of the log of the EMG profile
        EMG = ProcData.sleep.parameters.EMG{b};
        idx_raw = round(linspace(1,150,27));
        for ik = 1:1:26
%             if ik == 1
%             EMG_column(b,ik) = nanmean(EMG(1):EMG(idx_raw(ik)));
%             WhiskEvents_column(b,ik) = sum(WhiskEvents(1):WhiskEvents(idx_raw(ik)));
%             ForceEvents_column(b,ik) = sum(ForceEvents(1):ForceEvents(idx_raw(ik)));
%             else
            EMG_column(b,ik) = mean(EMG(idx_raw(ik):idx_raw(ik+1)));
            WhiskEvents_column(b,ik) = sum(WhiskEvents(idx_raw(ik):idx_raw(ik+1)));
            ForceEvents_column(b,ik) = sum(ForceEvents(idx_raw(ik):idx_raw(ik+1)));
%             end
        end
        
    end

   CortDeltaThetaLH = CortDeltaLH_column./CortThetaLH_column;
   CortDeltaThetaRH = CortDeltaRH_column./CortThetaRH_column;

   ModelDataFull = cell(length(CortDeltaLH_column),1);
    for kn = 1:1:length(CortDeltaLH_column)
    Dummy= [CortDeltaThetaLH(kn,:);CortDeltaThetaRH(kn,:);CortBetaLH_column(kn,:);CortBetaRH_column(kn,:);CortGammaLH_column(kn,:);CortGammaRH_column(kn,:);...
        HippDelta_column(kn,:); HippTheta_column(kn,:);HippGamma_column(kn,:);EMG_column(kn,:);WhiskEvents_column(kn,:);ForceEvents_column(kn,:)];
    ModelDataFull{kn} = Dummy';
    end
    ModelData.Features = ModelDataFull;
%         ModelData.Features = ModelDataFull;

    save(modelDataSetID,'ModelData')
end

end
