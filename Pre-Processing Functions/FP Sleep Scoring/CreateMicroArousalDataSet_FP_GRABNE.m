function [] = CreateMicroArousalDataSet_FP_GRABNE(procDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Md Shakhawat Hossin
% The Pennsylvania State University, Dept. of Biomedical Engineering
% Adopted from Kevin Turner
%________________________________________________________________________________________________________________________
%   Purpose: Go through each file and train a data set for the model or for model validation
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    MADataFileID = [procDataFileID(1:end-12) 'MAData.mat'];
    trainingDataFileID = [procDataFileID(1:end-12) 'TrainingData.mat'];
    if exist(trainingDataFileID,'file') 
        disp(['Loading ' procDataFileID ' for manual microarousal scoring.' ]); disp(' ')
        load(procDataFileID)

          if isfield(ProcData,'MicroArousals') == false ||  isfield(ProcData,'MicroArousals') == true
            load(trainingDataFileID)
            TableSize = height(trainingTable);
            DataSize = TableSize*5;
            MALabels = zeros(DataSize,1);
            
            OGLabels = trainingTable.behavState;
            for SL = 1:1:TableSize
                if OGLabels(SL) == "Not Sleep"
                    for ML = 1:1:5
                        MALabels(((SL-1)*5)+ML) = 1;
                    end
                elseif OGLabels(SL) == "NREM Sleep"
                    for ML = 1:1:5
                        MALabels(((SL-1)*5)+ML) = 2;
                    end
                elseif OGLabels(SL) == "REM Sleep"
                    for ML = 1:1:5
                        MALabels(((SL-1)*5)+ML) = 3;
                    end
                end
            end
            %% Binerized EMG
            EMG = abs(ProcData.data.EMG.emgSignal);
            ProcData.notes.EMGThresholdPerc = 95;
    
            EMGThreshold1 = prctile(EMG,ProcData.notes.EMGThresholdPerc)*ones(size(EMG));
            emgBinarized = EMG;
            emgBinarized(emgBinarized > EMGThreshold1) = 1;
            emgBinarized(emgBinarized < EMGThreshold1) = 0;
            emgBinarized(emgBinarized ==0 ) = nan;
    
            
            EMGCheck = figure; 
            subplot(211);plot((1:length(EMG))/ProcData.notes.dsFs,EMG,'color','k','LineWidth',1);      
            hold on; plot((1:length(EMG))/ProcData.notes.dsFs,-0.025*max(EMG)*emgBinarized,'|'); 
            xlim([0,ProcData.notes.trialDuration_sec])
            drawnow;
            hold off;
    
            ThresholdChoice = input('Do you want to change the threshold (Default is 95 percentile) for EMG scoring? ','s' ); disp(' ')
            if ThresholdChoice == 'y'
                change_Choice = 'y';
                while change_Choice == 'y'
                ProcData.notes.EMGThresholdPerc = input('Enter a threshold percentile: ');disp(' ')
                EMGThreshold1 = prctile(EMG,ProcData.notes.EMGThresholdPerc)*ones(size(EMG));
                emgBinarized = EMG;
                emgBinarized(emgBinarized > EMGThreshold1) = 1;
                emgBinarized(emgBinarized < EMGThreshold1) = 0;
                emgBinarized(emgBinarized ==0 ) = nan;
                subplot(212);plot((1:length(EMG))/ProcData.notes.dsFs,EMG,'color','k','LineWidth',1);      
                hold on; plot((1:length(EMG))/ProcData.notes.dsFs,-0.025*max(EMG)*emgBinarized,'|'); 
                xlim([0,ProcData.notes.trialDuration_sec])
                drawnow;
                hold off;
                change_Choice = input('Do you want to change the threshold? ','s');disp(' ')
                end
            end

            close(EMGCheck);
            
            EMGArousalLabels = zeros(DataSize,1);
            for EM = 1:1:DataSize
                if nansum(emgBinarized(1+((EM-1)*30):(EM*30))) > 2
                    EMGArousalLabels(EM) = 1;
                end
            end
    %         close(EMGCheck)     
            %% gamma power
            %{
            gammaPower = ProcData.data.cortical_LH.gammaBandPower;
            %derivative 
            gammaPower_Diff = diff(gammaPower);
            [z3,p3,k3] = butter(4,0.1/(ProcData.notes.dsFs/2),'low');
            [sos3,g3] = zp2sos(z3,p3,k3);
            filtgammaPower_Diff = filtfilt(sos3,g3,gammaPower_Diff);
            filtgammaPower_Diff(1:90) = filtgammaPower_Diff(91:180);
            % Threshold
            ProcData.notes.GammaThresholdPerc = 88;
            GammaThreshold = prctile(filtgammaPower_Diff,ProcData.notes.GammaThresholdPerc)*ones(size(filtgammaPower_Diff));
            GammaBinarized = filtgammaPower_Diff;
            GammaBinarized(GammaBinarized > GammaThreshold) = 1;
            GammaBinarized(GammaBinarized < GammaThreshold) = 0;
    %         GammaBinarized(GammaBinarized ==0 ) = nan;
            
            GammaArousalLabels = zeros(DataSize,1);
            for EM = 1:1:DataSize
                if EM == DataSize
                    if nansum(GammaBinarized(1+((EM-1)*30):(EM*30)-1)) > 3
                    GammaArousalLabels(EM) = 1;
                    end
                
                else
                    if nansum(GammaBinarized(1+((EM-1)*30):(EM*30))) > 3
                    GammaArousalLabels(EM) = 1;
                    end
            
                end
            
            end
            %}
            %% Detect microarousals and differentiate it from arousals
            % detect EMG triggered arousals.
            CMicroArousals = EMGArousalLabels ; %GammaArousalLabels | EMGArousalLabels ;
            MicroLabels =zeros(size(MALabels));
            for b = 1:DataSize
                if MALabels(b) == 2 % only detect during NREM Sleep
                    if CMicroArousals(b) == 1
                        MicroLabels(b) = 1;
                        continue;
                    end
                end
            end
            AwakeLabels = (MALabels == 1); % only take the not sleep labels from sleep score 
            MicroLabels = logical(MicroLabels); % convert to logical
            NMicroArousals = AwakeLabels | MicroLabels; % combined EMG deteced microarousals and user detected microarousals
            % matrix for better visualization
            NMatr(:,1) = double(NMicroArousals); % combined EMG deteced microarousals and user detected microarousals
            NMatr(:,2) = cumsum(NMicroArousals); % duration of any states
            NMatr_Arousal = [1; diff(NMicroArousals)]; % detect changes in states. either arousal or sleep
    
            ANT_A = (find(NMatr_Arousal == 1)); % positions where transition to arousals;
            ANT_S = (find(NMatr_Arousal == -1)); % positions where transition to arousals;
            
            if length(ANT_A) ~= length(ANT_S)
                ANT_A(end) = [];
            end
    
            A_N_T(:,1) = ANT_A; % positions where transition to arousals
            A_N_T(:,2) = ANT_S; % positions where transition from arousals
            A_N_T(:,3) = A_N_T(:,2) - A_N_T(:,1);   % duration of arousals states
           
            A_N_T(A_N_T(:,3) < 11 ,4)  = 1; % microarousals are arousals shorter than 10 seconds.

            MicroArousals_Final = zeros(size(AwakeLabels)); % final microarousals
            for AL = 1:1:length(A_N_T(:,3))
                if A_N_T(AL,4) == 1
                    MicroArousals_Final(A_N_T(AL,1):A_N_T(AL,2)-1) = 1; % create matrix for only microarousals
                end
            end
            NMatr(:,4) = MicroArousals_Final; % microarousals            
            %%
            saveFigs = 'y';
             [figHandle,ax1,ax2,ax3,ax4,ax5,ax7] = GenerateSingleFigures_MicroArousal_FP_GRABAchNE_(procDataFileID,saveFigs,MALabels,EMGArousalLabels,ProcData.notes,MicroArousals_Final);      
             % [figHandle,ax1,ax2,ax3] = GenerateSingleFigures_Proposal(procDataFileID,saveFigs,MALabels,EMGArousalLabels,ProcData.notes,MicroArousals_Final);      
    
%             [figHandle,ax1,ax2,ax3,ax4,ax5,ax6,ax7] = GenerateSingleFigures_MicroArousal_FP_GRABNE(procDataFileID,saveFigs,MALabels,EMGArousalLabels,ProcData.notes,GammaArousalLabels);      
            %%
            %{
            for b = 1:DataSize
                if MALabels(b) == 2
                    if CMicroArousals(b) == 1
                        MALabels(b) = 1;
                        MicroLabels(b) = 1;
                        continue;
                        
                    else % no EMG change detected. Need to look for microarousals without EMG activation
                        global buttonState %#ok<TLEV>
                                buttonState = 0;
                        global MArousal 
                               MArousal= 0; 
                        if b <= 3000
                        xaxis_1st = (ceil(b/300)-1)*300;
                        xaxis_2nd = (ceil(b/300))*300;
                        elseif b >= 3001
                        xaxis_1st = 3000;
                        xaxis_2nd = 3120;
                        end
    
                        xlim([xaxis_1st xaxis_2nd])
                
                        xStartVal = b;
                        xEndVal = b+1;
                        xInds = xStartVal:1:xEndVal;
                        figHandle = gcf;
            
                        subplot(ax3)
                        hold on
                        leftEdge3 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
                        rightEdge3 = xline(xInds(2),'color',colors('electric purple'),'LineWidth',2);
            
                        subplot(ax4)
                        hold on
                        leftEdge6 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
                        rightEdge6 = xline(xInds(2),'color',colors('electric purple'),'LineWidth',2);
        
                        MicroArousalSelection;
        
                        while buttonState == 0
                            drawnow()
                            if buttonState == 1
                                if MArousal == 1
                                    MALabels(b) = 1;
                                end
                                break;
                            end
                        end
        
                        delete(leftEdge3)
                        delete(leftEdge6)
                        delete(rightEdge3)
                        delete(rightEdge6)
                        
                    end
                end
            end
    %}
            close all
    %         ProcData.MicroArousal.GammaArousalLabels = GammaArousalLabels;
    %         ProcData.MicroArousal.EMGArousalLabels = EMGArousalLabels;
    %         ProcData.MicroArousal.CombinedMicroArousalLabels = CMicroArousals;
            ProcData.MicroArousals.MALogical = logical(MicroArousals_Final);
            save(MADataFileID, 'MALabels')
            save(procDataFileID, 'ProcData','-v7.3')
            clear NMatr A_N_T NMatr
          else
             disp([procDataFileID ' microArousal data exist. Continuing...']); disp(' ')
          end

    else
        disp([trainingDataFileID ' no training data exists.Please sleepscore the data first. Continuing...']); disp(' ')
    end
end

end
