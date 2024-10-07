function [] = CreateModelDataSet_FP_Shak(procDataFileIDs,NBins)
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
    CortDeltaLH_column = zeros(NBins,150);
    CortThetaLH_column = zeros(NBins,150);
    CortAlphaLH_column = zeros(NBins,150);
    CortBetaLH_column = zeros(NBins,150);
    CortGammaLH_column = zeros(NBins,150);

    CortDeltaRH_column = zeros(NBins,150);
    CortThetaRH_column = zeros(NBins,150);
    CortAlphaRH_column = zeros(NBins,150);
    CortBetaRH_column = zeros(NBins,150);
    CortGammaRH_column = zeros(NBins,150);

    HippDelta_column = zeros(NBins,150);
    HippTheta_column = zeros(NBins,150);
    HippGamma_column = zeros(NBins,150);
    WhiskEvents_column = zeros(NBins,150);
    ForceEvents_column = zeros(NBins,150);
    Whisk_column = zeros(NBins,150);
    Force_column = zeros(NBins,150);
    EMG_column = zeros(NBins,150);
%     avgHeartRate_column = zeros(NBins,1);
    % extract relevant parameters from each epoch
    for b = 1:length(CortDeltaLH_column)
        % Cortical delta
        CortDeltaLH_column(b,:) = ((ProcData.sleep.parameters.cortical_LH.deltaBandPower{b}));
        CortDeltaRH_column(b,:) = ((ProcData.sleep.parameters.cortical_RH.deltaBandPower{b}));
        % Cortical theta
        CortThetaLH_column(b,:) = ((ProcData.sleep.parameters.cortical_LH.thetaBandPower{b}));
        CortThetaRH_column(b,:) = ((ProcData.sleep.parameters.cortical_RH.thetaBandPower{b}));
        % Cortical alpha
        CortAlphaLH_column(b,:) = ((ProcData.sleep.parameters.cortical_LH.alphaBandPower{b}));
        CortAlphaRH_column(b,:) = ((ProcData.sleep.parameters.cortical_RH.alphaBandPower{b}));
        % Cortical beta
        CortBetaLH_column(b,:) = ((ProcData.sleep.parameters.cortical_LH.betaBandPower{b}));
        CortBetaRH_column(b,:) = ((ProcData.sleep.parameters.cortical_RH.betaBandPower{b}));
        % Cortical gamma
        CortGammaLH_column(b,:) = ((ProcData.sleep.parameters.cortical_LH.gammaBandPower{b}));
        CortGammaRH_column(b,:) = ((ProcData.sleep.parameters.cortical_RH.gammaBandPower{b}));
        % Hippocampal delta
        HippDelta_column(b,:) = ((ProcData.sleep.parameters.hippocampus.deltaBandPower{b}));
        % Hippocampal theta
        HippTheta_column(b,:) = ((ProcData.sleep.parameters.hippocampus.thetaBandPower{b}));
        % Hippocampal theta
        HippGamma_column(b,:) = ((ProcData.sleep.parameters.hippocampus.gammaBandPower{b}));
        if b==length(CortDeltaLH_column)
            % number of binarized whisking events
            WhiskEvents_column(b,:) = [(ProcData.sleep.parameters.binWhiskerAngle{b}),0,0];
            % number of binarized force sensor events
            ForceEvents_column(b,:) = [(ProcData.sleep.parameters.binForceSensor{b}),0];  
        else
            % number of binarized whisking events
            WhiskEvents_column(b,:) = [(ProcData.sleep.parameters.binWhiskerAngle{b})];
            % number of binarized force sensor events
            ForceEvents_column(b,:) = [(ProcData.sleep.parameters.binForceSensor{b})];
        end    
        % number of binarized whisking events
        Whisk_column(b,:) = [(ProcData.sleep.parameters.whiskerAngle{b})];
        % number of binarized force sensor events
        Force_column(b,:) = [(ProcData.sleep.parameters.ForceSensor{b})];
       
        % average of the log of the EMG profile
        EMG = ProcData.sleep.parameters.EMG{b};
        EMG_column(b,:) = (EMG);
    end
% 
%    CortDeltaThetaLH = CortDeltaLH_column./CortThetaLH_column;
%    CortDeltaThetaRH = CortDeltaRH_column./CortThetaRH_column;

   ModelDataFull = cell(length(CortDeltaLH_column),1);
    for kn = 1:1:length(CortDeltaLH_column)
%     Dummy= [CortDeltaThetaLH(kn,:);CortDeltaThetaRH(kn,:);CortBetaLH_column(kn,:);CortBetaRH_column(kn,:);CortGammaLH_column(kn,:);CortGammaRH_column(kn,:);...
%         HippDelta_column(kn,:); HippTheta_column(kn,:);HippGamma_column(kn,:);EMG_column(kn,:);WhiskEvents_column(kn,:);ForceEvents_column(kn,:)];
     Dummy = [CortDeltaLH_column(kn,:);CortDeltaRH_column(kn,:);CortThetaLH_column(kn,:);CortAlphaRH_column(kn,:);CortAlphaLH_column(kn,:);CortThetaRH_column(kn,:);CortBetaLH_column(kn,:);CortBetaRH_column(kn,:);CortGammaLH_column(kn,:);CortGammaRH_column(kn,:);...
        HippDelta_column(kn,:); HippTheta_column(kn,:);HippGamma_column(kn,:);...
        EMG_column(kn,:); Whisk_column(kn,:);Force_column(kn,:);...
        WhiskEvents_column(kn,:).*ForceEvents_column(kn,:).*EMG_column(kn,:);];
    
    ModelDataFull{kn} = Dummy';
    end
    ModelData.Features = ModelDataFull;
    save(modelDataSetID,'ModelData')
end

end
