function [AnalysisResults] = Fig1_S9_OptoStim_GRABNE_Response_CI_SleepDetected(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FPanimalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S2 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________
% load('S:\NEACh\AnalysisResults_firstHrs.mat');
% AnalysisResults = AnalysisResults_firstHrs;
FigInital = strfind(rootFolder,'\');
ManipulationType = rootFolder(FigInital(end)+1:end);
%% set-up and process data
solenoidNames = {'OptoStim'};
compDataTypes = {'Opto'};
OptoTypes = {'AllOptoStim' , 'NoSleep_OptoStim' , 'CauseArousal_OptoStim', 'CauseSleep_OptoStim'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
        for dd = 1:length(OptoTypes)
            OptoType = OptoTypes{1,dd};
            data.P_NE.OptoStim.(OptoType).count(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).P_NE.OptoStim.count;
            % mean
            data.P_NE.OptoStim.(OptoType).CBV(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).P_NE.OptoStim.CBV.CBV;
            data.P_ACh.OptoStim.(OptoType).CBV(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).P_ACh.OptoStim.CBV.CBV;
            data.P_NE.OptoStim.(OptoType).GFP(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).P_NE.OptoStim.GFP.GFP;
            data.P_ACh.OptoStim.(OptoType).GFP(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).P_ACh.OptoStim.GFP.GFP;
            %std
            data.P_NE.OptoStim.(OptoType).CBV_std(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).P_NE.OptoStim.CBV.CBVStD;
            data.P_ACh.OptoStim.(OptoType).CBV_std(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).P_ACh.OptoStim.CBV.CBVStD;
            data.P_NE.OptoStim.(OptoType).GFP_std(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).P_NE.OptoStim.GFP.GFPStD;
            data.P_ACh.OptoStim.(OptoType).GFP_std(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).P_ACh.OptoStim.GFP.GFPStD;
            %Pupil
            data.Pupil.OptoStim.(OptoType).Diameter(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).Pupil.OptoStim.Diameter.DiameterData;  
            data.Pupil.OptoStim.(OptoType).Diameter_std(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).Pupil.OptoStim.Diameter.DiameterStD;  
            % neural activity
            data.cortical_LH.OptoStim.(OptoType).cortMUA(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).cortical_LH.OptoStim.MUA.corticalData;
            data.cortical_LH.OptoStim.(OptoType).cortGam(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).cortical_LH.OptoStim.Gam.corticalData;
            data.cortical_LH.OptoStim.(OptoType).timeVector(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).cortical_LH.OptoStim.timeVector;
            data.cortical_LH.OptoStim.(OptoType).cortS_LH(:,:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).cortical_LH.OptoStim.ECoG.corticalS_LH;
            data.cortical_LH.OptoStim.(OptoType).cortS_LH_Gam(:,:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).cortical_LH.OptoStim.ECoG.corticalS_LH(49:end,20:23);
            % data.cortical_RH.OptoStim.(OptoType).cortS_RH(:,:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).cortical_LH.OptoStim.ECoG.corticalS_RH;
            % data.cortical_RH.OptoStim.(OptoType).cortS_RH_Gam(:,:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).cortical_LH.OptoStim.ECoG.corticalS_RH(49:end,20:23);
            data.cortical_LH.OptoStim.(OptoType).T(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).cortical_LH.OptoStim.ECoG.T;
            data.cortical_LH.OptoStim.(OptoType).F(:,aa) = AnalysisResults.(animalID).OptoStim.(OptoType).cortical_LH.OptoStim.ECoG.F;
        end
end
% extract the raw signals
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
        for dd = 1:length(OptoTypes)
            OptoType = OptoTypes{1,dd};
            % the size of the matrix
            MatLength = size(AnalysisResults.(animalID).OptoStim.(OptoType).P_NE.OptoStim.CBV.CBVRaw,1);
            if aa == 1
                DataLength = 0;
            else
                DataLength = size(data.P_NE.OptoStim.(OptoType).CBVRaw,1);
            end
            % RawData
            data.P_NE.OptoStim.(OptoType).CBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).OptoStim.(OptoType).P_NE.OptoStim.CBV.CBVRaw;
            data.P_ACh.OptoStim.(OptoType).CBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).OptoStim.(OptoType).P_ACh.OptoStim.CBV.CBVRaw;
            data.P_NE.OptoStim.(OptoType).GFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).OptoStim.(OptoType).P_NE.OptoStim.GFP.GFPRaw;
            data.P_ACh.OptoStim.(OptoType).GFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).OptoStim.(OptoType).P_ACh.OptoStim.GFP.GFPRaw;
            data.Pupil.OptoStim.(OptoType).DiameterRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).OptoStim.(OptoType).Pupil.OptoStim.Diameter.DiameterRaw;
        end
end

%% concatenate the data

for ff = 1:length(OptoTypes)
    OptoType = OptoTypes{1,ff};
    data.(OptoType).count = data.P_NE.OptoStim.(OptoType).count;
    
    data.(OptoType).P_AChCBV = data.P_ACh.OptoStim.(OptoType).CBV;
    data.(OptoType).P_AChGFP = data.P_ACh.OptoStim.(OptoType).GFP;
    data.(OptoType).P_NECBV = data.P_NE.OptoStim.(OptoType).CBV;
    data.(OptoType).P_NEGFP = data.P_NE.OptoStim.(OptoType).GFP;
    
    data.(OptoType).P_AChCBV_std = data.P_ACh.OptoStim.(OptoType).CBV_std;
    data.(OptoType).P_AChGFP_std = data.P_ACh.OptoStim.(OptoType).GFP_std;
    data.(OptoType).P_NECBV_std = data.P_NE.OptoStim.(OptoType).CBV_std;
    data.(OptoType).P_NEGFP_std = data.P_NE.OptoStim.(OptoType).GFP_std;
    
    data.(OptoType).PupilDiameter = data.Pupil.OptoStim.(OptoType).Diameter;
    data.(OptoType).PupilDiameter_std = data.Pupil.OptoStim.(OptoType).Diameter_std;
    
    data.(OptoType).Cortical.cortMUA = data.cortical_LH.OptoStim.(OptoType).cortMUA; 
    data.(OptoType).Cortical.cortGam = data.cortical_LH.OptoStim.(OptoType).cortGam; 
    data.(OptoType).Cortical.timeVector = data.cortical_LH.OptoStim.(OptoType).timeVector;
    data.(OptoType).Cortical.cortS_LH = data.cortical_LH.OptoStim.(OptoType).cortS_LH; 
    data.(OptoType).Cortical.cortS_LH_Gam = data.cortical_LH.OptoStim.(OptoType).cortS_LH_Gam;
    % data.(OptoType).Cortical.cortS_RH = data.cortical_LH.OptoStim.(OptoType).cortS_RH; 
    % data.(OptoType).Cortical.cortS_RH_Gam = data.cortical_LH.OptoStim.(OptoType).cortS_RH_Gam;
    data.(OptoType).Cortical.T = data.cortical_LH.OptoStim.(OptoType).T;
    data.(OptoType).Cortical.F = data.cortical_LH.OptoStim.(OptoType).F;
end

%% take the averages of each field through the proper dimension
for ff = 1:length(OptoTypes)
    OptoType = OptoTypes{1,ff};
    data.(OptoType).mean_Count = mean(data.(OptoType).count,2);
    data.(OptoType).std_Count = std(data.(OptoType).count,0,2);   
    
    data.(OptoType).P_AChmean_CBV = mean(data.(OptoType).P_AChCBV,2);
    data.(OptoType).P_AChstd_CBV = std(data.(OptoType).P_AChCBV,0,2);

    data.(OptoType).P_NEmean_CBV = mean(data.(OptoType).P_NECBV,2);
    data.(OptoType).P_NEstd_CBV = std(data.(OptoType).P_NECBV,0,2);

    data.(OptoType).P_AChmean_GFP = mean(data.(OptoType).P_AChGFP,2);
    data.(OptoType).P_AChstd_GFP = std(data.(OptoType).P_AChGFP,0,2);

    data.(OptoType).P_NEmean_GFP = mean(data.(OptoType).P_NEGFP,2);
    data.(OptoType).P_NEstd_GFP = std(data.(OptoType).P_NEGFP,0,2);

    data.(OptoType).P_AChmean_CBV_std = mean(data.(OptoType).P_AChCBV_std,2);
    data.(OptoType).P_NEmean_CBV_std = mean(data.(OptoType).P_NECBV_std,2);

    data.(OptoType).P_AChmean_GFP_std = mean(data.(OptoType).P_AChGFP_std,2);
    data.(OptoType).P_NEmean_GFP_std = mean(data.(OptoType).P_NEGFP_std,2);

    data.(OptoType).Pupilmean_Diameter = mean(data.(OptoType).PupilDiameter,2);
    data.(OptoType).Pupilstd_Diameter = std(data.(OptoType).PupilDiameter,0,2);

    data.(OptoType).mean_CortMUA = mean(data.(OptoType).Cortical.cortMUA,2);
    data.(OptoType).std_CortMUA = std(data.(OptoType).Cortical.cortMUA,0,2);
    data.(OptoType).mean_CortGam = mean(data.(OptoType).Cortical.cortGam,2);
    data.(OptoType).std_CortGam = std(data.(OptoType).Cortical.cortGam,0,2);
    data.(OptoType).mean_timeVector = mean(data.(OptoType).Cortical.timeVector,2);

    data.(OptoType).mean_CortS_LH = mean(data.(OptoType).Cortical.cortS_LH,3).*100;
    data.(OptoType).mean_CortS_LH_Gam = mean(mean(mean(data.(OptoType).Cortical.cortS_LH_Gam.*100,2),2),3);
    data.(OptoType).std_CortS_LH_Gam = std(mean(mean(data.(OptoType).Cortical.cortS_LH_Gam.*100,2),2),0,3);

    % data.(OptoType).mean_CortS_RH = mean(data.(OptoType).Cortical.cortS_RH,3).*100;
    % data.(OptoType).mean_CortS_RH_Gam = mean(mean(mean(data.(OptoType).Cortical.cortS_RH_Gam.*100,2),2),3);
    % data.(OptoType).std_CortS_RH_Gam = std(mean(mean(data.(OptoType).Cortical.cortS_RH_Gam.*100,2),2),0,3);

    data.(OptoType).mean_T = mean(data.(OptoType).Cortical.T,2);
    data.(OptoType).mean_F = mean(data.(OptoType).Cortical.F,2);
end

%% raw data
for ff = 1:length(OptoTypes)
    OptoType = OptoTypes{1,ff};
    data.(OptoType).P_AChCBVRaw = data.P_ACh.OptoStim.(OptoType).CBVRaw;
    data.(OptoType).P_AChGFPRaw = data.P_ACh.OptoStim.(OptoType).GFPRaw;
    data.(OptoType).P_NECBVRaw = data.P_NE.OptoStim.(OptoType).CBVRaw;
    data.(OptoType).P_NEGFPRaw = data.P_NE.OptoStim.(OptoType).GFPRaw;
    data.(OptoType).PupilDiameterRaw = data.Pupil.OptoStim.(OptoType).DiameterRaw;
end

% take the averages of each field through the proper dimension
for ff = 1:length(OptoTypes)
    OptoType = OptoTypes{1,ff};
    data.(OptoType).P_AChCBV_N_Exp = size(data.(OptoType).P_AChCBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(OptoType).P_AChCBV_Mean = mean(data.(OptoType).P_AChCBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(OptoType).P_AChCBV_SEM = std(data.(OptoType).P_AChCBVRaw,1)/sqrt(data.(OptoType).P_AChCBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(OptoType).P_AChCBV_CI95 = tinv([0.025 0.975], data.(OptoType).P_AChCBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(OptoType).P_AChCBV_yCI95 = bsxfun(@times, data.(OptoType).P_AChCBV_SEM, data.(OptoType).P_AChCBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(OptoType).P_NECBV_N_Exp = size(data.(OptoType).P_NECBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(OptoType).P_NECBV_Mean = mean(data.(OptoType).P_NECBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(OptoType).P_NECBV_SEM = std(data.(OptoType).P_NECBVRaw,1)/sqrt(data.(OptoType).P_NECBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(OptoType).P_NECBV_CI95 = tinv([0.025 0.975], data.(OptoType).P_NECBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(OptoType).P_NECBV_yCI95 = bsxfun(@times, data.(OptoType).P_NECBV_SEM, data.(OptoType).P_NECBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(OptoType).P_AChGFP_N_Exp = size(data.(OptoType).P_AChGFPRaw,1);
    data.(OptoType).P_AChGFP_Mean = mean(data.(OptoType).P_AChGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(OptoType).P_AChGFP_SEM = std(data.(OptoType).P_AChGFPRaw,1)/sqrt(data.(OptoType).P_AChGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(OptoType).P_AChGFP_CI95 = tinv([0.025 0.975], data.(OptoType).P_AChGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(OptoType).P_AChGFP_yCI95 = bsxfun(@times, data.(OptoType).P_AChGFP_SEM, data.(OptoType).P_AChGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(OptoType).P_NEGFP_N_Exp = size(data.(OptoType).P_NEGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(OptoType).P_NEGFP_Mean = mean(data.(OptoType).P_NEGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(OptoType).P_NEGFP_SEM = std(data.(OptoType).P_NEGFPRaw,1)/sqrt(data.(OptoType).P_NEGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(OptoType).P_NEGFP_CI95 = tinv([0.025 0.975], data.(OptoType).P_NEGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(OptoType).P_NEGFP_yCI95 = bsxfun(@times, data.(OptoType).P_NEGFP_SEM, data.(OptoType).P_NEGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(OptoType).PupilDiameter_N_Exp = size(data.(OptoType).PupilDiameterRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(OptoType).PupilDiameter_Mean = mean(data.(OptoType).PupilDiameterRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(OptoType).PupilDiameter_SEM = std(data.(OptoType).PupilDiameterRaw,1)/sqrt(data.(OptoType).PupilDiameter_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(OptoType).PupilDiameter_CI95 = tinv([0.025 0.975], data.(OptoType).PupilDiameter_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(OptoType).PupilDiameter_yCI95 = bsxfun(@times, data.(OptoType).PupilDiameter_SEM, data.(OptoType).PupilDiameter_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    % data.(OptoType).P_CBV_N_Exp = size(data.(OptoType).P_CBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    % data.(OptoType).P_CBV_Mean = mean(data.(OptoType).P_CBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    % data.(OptoType).P_CBV_SEM = std(data.(OptoType).P_CBVRaw,1)/sqrt(data.(OptoType).P_CBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    % data.(OptoType).P_CBV_CI95 = tinv([0.025 0.975], data.(OptoType).P_CBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    % data.(OptoType).P_CBV_yCI95 = bsxfun(@times, data.(OptoType).P_CBV_SEM, data.(OptoType).P_CBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

end
%% plot the response
summaryFigureN = figure ;
%% All Opto
ax1 = subplot(2,2,1);
p1= plot(data.AllOptoStim.mean_timeVector,data.AllOptoStim.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on
NEn =  data.AllOptoStim.P_AChCBV_Mean + data.AllOptoStim.P_AChCBV_yCI95;
patch([data.AllOptoStim.mean_timeVector' fliplr(data.AllOptoStim.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0.8500 0.3250 0.0980],'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.AllOptoStim.mean_timeVector,data.AllOptoStim.P_AChCBV_Mean + data.AllOptoStim.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)

p2= plot(data.AllOptoStim.mean_timeVector,data.AllOptoStim.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
NEn =  data.AllOptoStim.P_NECBV_Mean + data.AllOptoStim.P_NECBV_yCI95;
patch([data.AllOptoStim.mean_timeVector' fliplr(data.AllOptoStim.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.AllOptoStim.mean_timeVector,data.AllOptoStim.P_NECBV_Mean + data.AllOptoStim.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
ylabel('CBV \DeltaF/F (%)')
ylim([-1.5 1]);

yyaxis right
p3= plot(data.AllOptoStim.mean_timeVector,data.AllOptoStim.PupilDiameter_Mean,'-','color','k','LineWidth',2);
hold on
NEn =  data.AllOptoStim.PupilDiameter_Mean + data.AllOptoStim.PupilDiameter_yCI95;
patch([data.AllOptoStim.mean_timeVector' fliplr(data.AllOptoStim.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],'k','FaceColor','k','EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.AllOptoStim.mean_timeVector,data.AllOptoStim.PupilDiameter_Mean + data.AllOptoStim.PupilDiameter_yCI95,'-','color','k','LineWidth',0.10)
ylabel('Diameter (Z)')

xline(0,'-',{'Start'});
xline(5,'-',{'End'});
title('All optostims')
ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = 'k';
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'LH', 'RH','Pupil','Location','best')
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ylim([-0.5 1]);
xlim([-5 15])
axis square

% ECoG Spectrogram LH
% ax2 = subplot(4,2,3);
% imagesc(data.AllOptoStim.mean_T,data.AllOptoStim.mean_F,data.AllOptoStim.mean_CortS_LH)
% title('All optostims')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c4 = colorbar;
% ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% % yscale log
% axis square
% axis xy
% set(gca,'box','off')
% % xlim([-5 15])
% xline(0,'-',{'Start'});
% xline(3,'-',{'End'});

%% No Sleep Opto
ax3 = subplot(2,2,2);
p1= plot(data.NoSleep_OptoStim.mean_timeVector,data.NoSleep_OptoStim.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on
NEn =  data.NoSleep_OptoStim.P_AChCBV_Mean + data.NoSleep_OptoStim.P_AChCBV_yCI95;
patch([data.NoSleep_OptoStim.mean_timeVector' fliplr(data.NoSleep_OptoStim.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0.8500 0.3250 0.0980],'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.NoSleep_OptoStim.mean_timeVector,data.NoSleep_OptoStim.P_AChCBV_Mean + data.NoSleep_OptoStim.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)

p2= plot(data.NoSleep_OptoStim.mean_timeVector,data.NoSleep_OptoStim.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
NEn =  data.NoSleep_OptoStim.P_NECBV_Mean + data.NoSleep_OptoStim.P_NECBV_yCI95;
patch([data.NoSleep_OptoStim.mean_timeVector' fliplr(data.NoSleep_OptoStim.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.NoSleep_OptoStim.mean_timeVector,data.NoSleep_OptoStim.P_NECBV_Mean + data.NoSleep_OptoStim.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
ylabel('CBV \DeltaF/F (%)')
ylim([-1.5 1]);

yyaxis right
p3= plot(data.NoSleep_OptoStim.mean_timeVector,data.NoSleep_OptoStim.PupilDiameter_Mean,'-','color','k','LineWidth',2);
hold on
NEn =  data.NoSleep_OptoStim.PupilDiameter_Mean + data.NoSleep_OptoStim.PupilDiameter_yCI95;
patch([data.NoSleep_OptoStim.mean_timeVector' fliplr(data.NoSleep_OptoStim.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],'k','FaceColor','k','EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.NoSleep_OptoStim.mean_timeVector,data.NoSleep_OptoStim.PupilDiameter_Mean + data.NoSleep_OptoStim.PupilDiameter_yCI95,'-','color','k','LineWidth',0.10)
ylabel('Diameter (Z)')

xline(0,'-',{'Start'});
xline(5,'-',{'End'});
title('Optostim During Arousal')

ax3.YAxis(1).Color = 'k';
ax3.YAxis(2).Color = 'k';
xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'LH', 'RH','Pupil','Location','best')
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
ylim([-0.5 1]);
xlim([-5 15])
axis square

% ECoG Spectrogram LH
% ax4 = subplot(4,2,4);
% imagesc(data.NoSleep_OptoStim.mean_T,data.NoSleep_OptoStim.mean_F,data.NoSleep_OptoStim.mean_CortS_LH)
% title('Optostim during Awake')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c4 = colorbar;
% ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% % yscale log
% axis square
% axis xy
% set(gca,'box','off')
% % xlim([-5 15])
% xline(0,'-',{'Start'});
% xline(3,'-',{'End'});

%% Opto Cause Arousal
ax5 = subplot(2,2,3);
p1= plot(data.CauseArousal_OptoStim.mean_timeVector,data.CauseArousal_OptoStim.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on
NEn =  data.CauseArousal_OptoStim.P_AChCBV_Mean + data.CauseArousal_OptoStim.P_AChCBV_yCI95;
patch([data.CauseArousal_OptoStim.mean_timeVector' fliplr(data.CauseArousal_OptoStim.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0.8500 0.3250 0.0980],'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.CauseArousal_OptoStim.mean_timeVector,data.CauseArousal_OptoStim.P_AChCBV_Mean + data.CauseArousal_OptoStim.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)

p2= plot(data.CauseArousal_OptoStim.mean_timeVector,data.CauseArousal_OptoStim.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
NEn =  data.CauseArousal_OptoStim.P_NECBV_Mean + data.CauseArousal_OptoStim.P_NECBV_yCI95;
patch([data.CauseArousal_OptoStim.mean_timeVector' fliplr(data.CauseArousal_OptoStim.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.CauseArousal_OptoStim.mean_timeVector,data.CauseArousal_OptoStim.P_NECBV_Mean + data.CauseArousal_OptoStim.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
ylabel('CBV \DeltaF/F (%)')
ylim([-1.5 1]);

yyaxis right
p3= plot(data.CauseArousal_OptoStim.mean_timeVector,data.CauseArousal_OptoStim.PupilDiameter_Mean,'-','color','k','LineWidth',2);
hold on
NEn =  data.CauseArousal_OptoStim.PupilDiameter_Mean + data.CauseArousal_OptoStim.PupilDiameter_yCI95;
patch([data.CauseArousal_OptoStim.mean_timeVector' fliplr(data.CauseArousal_OptoStim.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],'k','FaceColor','k','EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.CauseArousal_OptoStim.mean_timeVector,data.CauseArousal_OptoStim.PupilDiameter_Mean + data.CauseArousal_OptoStim.PupilDiameter_yCI95,'-','color','k','LineWidth',0.10)
ylabel('Diameter (Z)')

xline(0,'-',{'Start'});
xline(5,'-',{'End'});
title('Optostim Cause Arousal')

ax5.YAxis(1).Color = 'k';
ax5.YAxis(2).Color = 'k';
xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'LH', 'RH','Pupil','Location','best')
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
ylim([-0.5 1]);
xlim([-5 15])
axis square

% ECoG Spectrogram LH
% ax6 = subplot(4,2,7);
% imagesc(data.CauseArousal_OptoStim.mean_T,data.CauseArousal_OptoStim.mean_F,data.CauseArousal_OptoStim.mean_CortS_LH)
% title('Optostim Cause Arousal')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c4 = colorbar;
% ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% % yscale log
% axis square
% axis xy
% set(gca,'box','off')
% % xlim([-5 15])
% xline(0,'-',{'Start'});
% xline(3,'-',{'End'});
%% Opto Does not Cause Arousal
ax7 = subplot(2,2,4);
p1= plot(data.CauseSleep_OptoStim.mean_timeVector,data.CauseSleep_OptoStim.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on
NEn =  data.CauseSleep_OptoStim.P_AChCBV_Mean + data.CauseSleep_OptoStim.P_AChCBV_yCI95;
patch([data.CauseSleep_OptoStim.mean_timeVector' fliplr(data.CauseSleep_OptoStim.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0.8500 0.3250 0.0980],'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.CauseSleep_OptoStim.mean_timeVector,data.CauseSleep_OptoStim.P_AChCBV_Mean + data.CauseSleep_OptoStim.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)

p2= plot(data.CauseSleep_OptoStim.mean_timeVector,data.CauseSleep_OptoStim.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
NEn =  data.CauseSleep_OptoStim.P_NECBV_Mean + data.CauseSleep_OptoStim.P_NECBV_yCI95;
patch([data.CauseSleep_OptoStim.mean_timeVector' fliplr(data.CauseSleep_OptoStim.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.CauseSleep_OptoStim.mean_timeVector,data.CauseSleep_OptoStim.P_NECBV_Mean + data.CauseSleep_OptoStim.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
ylabel('CBV \DeltaF/F (%)')
ylim([-1.5 1]);

yyaxis right
p3= plot(data.CauseSleep_OptoStim.mean_timeVector,data.CauseSleep_OptoStim.PupilDiameter_Mean,'-','color','k','LineWidth',2);
hold on
NEn =  data.CauseSleep_OptoStim.PupilDiameter_Mean + data.CauseSleep_OptoStim.PupilDiameter_yCI95;
patch([data.CauseSleep_OptoStim.mean_timeVector' fliplr(data.CauseSleep_OptoStim.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],'k','FaceColor','k','EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.CauseSleep_OptoStim.mean_timeVector,data.CauseSleep_OptoStim.PupilDiameter_Mean + data.CauseSleep_OptoStim.PupilDiameter_yCI95,'-','color','k','LineWidth',0.10)
ylabel('Diameter (Z)')

xline(0,'-',{'Start'});
xline(5,'-',{'End'});
title('Optostim No Arousal')

ax7.YAxis(1).Color = 'k';
ax7.YAxis(2).Color = 'k';
xlabel('Peri-stimulus time (s)')
% legend([p1,p2,p3],'LH', 'RH','Pupil','Location','best')
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
ylim([-0.5 1]);
xlim([-5 15])
axis square

% ECoG Spectrogram LH
% ax8 = subplot(4,2,8);
% imagesc(data.CauseSleep_OptoStim.mean_T,data.CauseSleep_OptoStim.mean_F,data.CauseSleep_OptoStim.mean_CortS_LH)
% title('Optostim No Arousal')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c4 = colorbar;
% ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% % yscale log
% axis square
% axis xy
% set(gca,'box','off')
% % xlim([-5 15])
% xline(0,'-',{'Start'});
% xline(3,'-',{'End'});
%% Axes properties
% linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8],'x')
linkaxes([ax1,ax3,ax5,ax7],'x')

ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
% ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
% ax8Pos = get(ax8,'position');

% ax2Pos(3:4) = ax1Pos(3:4);
% ax4Pos(3:4) = ax3Pos(3:4);
% ax6Pos(3:4) = ax1Pos(3:4);
% ax8Pos(3:4) = ax3Pos(3:4);


set(ax1,'position',ax1Pos);
% set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
% set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
% set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
% set(ax8,'position',ax8Pos);
    %% save figure(s)
    if strcmp(saveFigs,'y') == true
        if firstHrs == "false"
            dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
        elseif firstHrs == "true"
            dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
        end    
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(summaryFigureN,[dirpath ManipulationType '_' 'OptoStim_Response_CI_ACh_NE']); %animalID   
        set(summaryFigureN,'PaperPositionMode','auto');
        print('-painters','-dpdf','-fillpage',[dirpath ManipulationType '_'  'OptoStim_Response_CI_ACh_NE']) % animalID _NoSleep NoSleep_
        close 
    end
end
