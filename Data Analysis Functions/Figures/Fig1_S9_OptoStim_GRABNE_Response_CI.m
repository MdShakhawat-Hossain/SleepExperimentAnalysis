function [AnalysisResults] = Fig1_S9_OptoStim_GRABNE_Response_CI(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FPanimalIDs)
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
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
        for dd = 1:length(solenoidNames)
            solenoidName = solenoidNames{1,dd};
            data.P_NE.(solenoidName).count(:,aa) = AnalysisResults.(animalID).OptoStim.P_NE.(solenoidName).count;
            % mean
            data.P_NE.(solenoidName).CBV(:,aa) = AnalysisResults.(animalID).OptoStim.P_NE.(solenoidName).CBV.CBV;
            data.P_ACh.(solenoidName).CBV(:,aa) = AnalysisResults.(animalID).OptoStim.P_ACh.(solenoidName).CBV.CBV;
            data.P_NE.(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).OptoStim.P_NE.(solenoidName).GFP.GFP;
            data.P_ACh.(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).OptoStim.P_ACh.(solenoidName).GFP.GFP;
            %std
            data.P_NE.(solenoidName).CBV_std(:,aa) = AnalysisResults.(animalID).OptoStim.P_NE.(solenoidName).CBV.CBVStD;
            data.P_ACh.(solenoidName).CBV_std(:,aa) = AnalysisResults.(animalID).OptoStim.P_ACh.(solenoidName).CBV.CBVStD;
            data.P_NE.(solenoidName).GFP_std(:,aa) = AnalysisResults.(animalID).OptoStim.P_NE.(solenoidName).GFP.GFPStD;
            data.P_ACh.(solenoidName).GFP_std(:,aa) = AnalysisResults.(animalID).OptoStim.P_ACh.(solenoidName).GFP.GFPStD;
            %Pupil
            data.Pupil.(solenoidName).Diameter(:,aa) = AnalysisResults.(animalID).OptoStim.Pupil.(solenoidName).Diameter.DiameterData;  
            data.Pupil.(solenoidName).Diameter_std(:,aa) = AnalysisResults.(animalID).OptoStim.Pupil.(solenoidName).Diameter.DiameterStD;  
            % neural activity
            data.cortical_LH.(solenoidName).cortMUA(:,aa) = AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoidName).MUA.corticalData;
            data.cortical_LH.(solenoidName).cortGam(:,aa) = AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoidName).Gam.corticalData;
            data.cortical_LH.(solenoidName).timeVector(:,aa) = AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoidName).timeVector;
            data.cortical_LH.(solenoidName).cortS_LH(:,:,aa) = AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoidName).ECoG.corticalS_LH;
            data.cortical_LH.(solenoidName).cortS_LH_Gam(:,:,aa) = AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoidName).ECoG.corticalS_LH(49:end,20:23);
            % data.cortical_RH.(solenoidName).cortS_RH(:,:,aa) = AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoidName).ECoG.corticalS_RH;
            % data.cortical_RH.(solenoidName).cortS_RH_Gam(:,:,aa) = AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoidName).ECoG.corticalS_RH(49:end,20:23);
            data.cortical_LH.(solenoidName).T(:,aa) = AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoidName).ECoG.T;
            data.cortical_LH.(solenoidName).F(:,aa) = AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoidName).ECoG.F;
        end
end
% extract the raw signals
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
        for dd = 1:length(solenoidNames)
            solenoidName = solenoidNames{1,dd};
            % the size of the matrix
            MatLength = size(AnalysisResults.(animalID).OptoStim.P_NE.(solenoidName).CBV.CBVRaw,1);
            if aa == 1
                DataLength = 0;
            else
                DataLength = size(data.P_NE.(solenoidName).CBVRaw,1);
            end
            % mean
            data.P_NE.(solenoidName).CBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).OptoStim.P_NE.(solenoidName).CBV.CBVRaw;
            data.P_ACh.(solenoidName).CBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).OptoStim.P_ACh.(solenoidName).CBV.CBVRaw;
            data.P_NE.(solenoidName).GFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).OptoStim.P_NE.(solenoidName).GFP.GFPRaw;
            data.P_ACh.(solenoidName).GFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).OptoStim.P_ACh.(solenoidName).GFP.GFPRaw;
            data.Pupil.(solenoidName).DiameterRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).OptoStim.Pupil.(solenoidName).Diameter.DiameterRaw;
        end
end

% concatenate the data from the contra and ipsi data
% contra
data.Opto.count = data.P_NE.OptoStim.count;
data.Opto.P_AChCBV = data.P_ACh.OptoStim.CBV;
data.Opto.P_AChGFP = data.P_ACh.OptoStim.GFP;
data.Opto.P_NECBV = data.P_NE.OptoStim.CBV;
data.Opto.P_NEGFP = data.P_NE.OptoStim.GFP;

data.Opto.P_AChCBV_std = data.P_ACh.OptoStim.CBV_std;
data.Opto.P_AChGFP_std = data.P_ACh.OptoStim.GFP_std;
data.Opto.P_NECBV_std = data.P_NE.OptoStim.CBV_std;
data.Opto.P_NEGFP_std = data.P_NE.OptoStim.GFP_std;

data.Opto.PupilDiameter = data.Pupil.OptoStim.Diameter;
data.Opto.PupilDiameter_std = data.Pupil.OptoStim.Diameter_std;

data.Opto.Cortical.cortMUA = data.cortical_LH.OptoStim.cortMUA; 
data.Opto.Cortical.cortGam = data.cortical_LH.OptoStim.cortGam; 
data.Opto.Cortical.timeVector = data.cortical_LH.OptoStim.timeVector;
data.Opto.Cortical.cortS_LH = data.cortical_LH.OptoStim.cortS_LH; 
data.Opto.Cortical.cortS_LH_Gam = data.cortical_LH.OptoStim.cortS_LH_Gam;
% data.Opto.Cortical.cortS_RH = data.cortical_LH.OptoStim.cortS_RH; 
% data.Opto.Cortical.cortS_RH_Gam = data.cortical_LH.OptoStim.cortS_RH_Gam;
data.Opto.Cortical.T = data.cortical_LH.OptoStim.T;
data.Opto.Cortical.F = data.cortical_LH.OptoStim.F;

% take the averages of each field through the proper dimension
for ff = 1:length(compDataTypes)
    compDataType = compDataTypes{1,ff};
    data.(compDataType).mean_Count = mean(data.(compDataType).count,2);
    data.(compDataType).std_Count = std(data.(compDataType).count,0,2);   
    
    data.(compDataType).P_AChmean_CBV = mean(data.(compDataType).P_AChCBV,2);
    data.(compDataType).P_AChstd_CBV = std(data.(compDataType).P_AChCBV,0,2);

    data.(compDataType).P_NEmean_CBV = mean(data.(compDataType).P_NECBV,2);
    data.(compDataType).P_NEstd_CBV = std(data.(compDataType).P_NECBV,0,2);

    data.(compDataType).P_AChmean_GFP = mean(data.(compDataType).P_AChGFP,2);
    data.(compDataType).P_AChstd_GFP = std(data.(compDataType).P_AChGFP,0,2);

    data.(compDataType).P_NEmean_GFP = mean(data.(compDataType).P_NEGFP,2);
    data.(compDataType).P_NEstd_GFP = std(data.(compDataType).P_NEGFP,0,2);

    data.(compDataType).P_AChmean_CBV_std = mean(data.(compDataType).P_AChCBV_std,2);
    data.(compDataType).P_NEmean_CBV_std = mean(data.(compDataType).P_NECBV_std,2);

    data.(compDataType).P_AChmean_GFP_std = mean(data.(compDataType).P_AChGFP_std,2);
    data.(compDataType).P_NEmean_GFP_std = mean(data.(compDataType).P_NEGFP_std,2);

    data.(compDataType).Pupilmean_Diameter = mean(data.(compDataType).PupilDiameter,2);
    data.(compDataType).Pupilstd_Diameter = std(data.(compDataType).PupilDiameter,0,2);

    data.(compDataType).mean_CortMUA = mean(data.(compDataType).Cortical.cortMUA,2);
    data.(compDataType).std_CortMUA = std(data.(compDataType).Cortical.cortMUA,0,2);
    data.(compDataType).mean_CortGam = mean(data.(compDataType).Cortical.cortGam,2);
    data.(compDataType).std_CortGam = std(data.(compDataType).Cortical.cortGam,0,2);
    data.(compDataType).mean_timeVector = mean(data.(compDataType).Cortical.timeVector,2);

    data.(compDataType).mean_CortS_LH = mean(data.(compDataType).Cortical.cortS_LH,3).*100;
    data.(compDataType).mean_CortS_LH_Gam = mean(mean(mean(data.(compDataType).Cortical.cortS_LH_Gam.*100,2),2),3);
    data.(compDataType).std_CortS_LH_Gam = std(mean(mean(data.(compDataType).Cortical.cortS_LH_Gam.*100,2),2),0,3);

    % data.(compDataType).mean_CortS_RH = mean(data.(compDataType).Cortical.cortS_RH,3).*100;
    % data.(compDataType).mean_CortS_RH_Gam = mean(mean(mean(data.(compDataType).Cortical.cortS_RH_Gam.*100,2),2),3);
    % data.(compDataType).std_CortS_RH_Gam = std(mean(mean(data.(compDataType).Cortical.cortS_RH_Gam.*100,2),2),0,3);

    data.(compDataType).mean_T = mean(data.(compDataType).Cortical.T,2);
    data.(compDataType).mean_F = mean(data.(compDataType).Cortical.F,2);
end
%% raw data
% concatenate the data from the contra and ipsi data
% contra
data.Opto.P_AChCBVRaw = data.P_ACh.OptoStim.CBVRaw;
data.Opto.P_AChGFPRaw = data.P_ACh.OptoStim.GFPRaw;
data.Opto.P_NECBVRaw = data.P_NE.OptoStim.CBVRaw;
data.Opto.P_NEGFPRaw = data.P_NE.OptoStim.GFPRaw;
data.Opto.PupilDiameterRaw = data.Pupil.OptoStim.DiameterRaw;

% take the averages of each field through the proper dimension
for ff = 1:length(compDataTypes)
    compDataType = compDataTypes{1,ff};
    data.(compDataType).P_AChCBV_N_Exp = size(data.(compDataType).P_AChCBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(compDataType).P_AChCBV_Mean = mean(data.(compDataType).P_AChCBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(compDataType).P_AChCBV_SEM = std(data.(compDataType).P_AChCBVRaw,1)/sqrt(data.(compDataType).P_AChCBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(compDataType).P_AChCBV_CI95 = tinv([0.025 0.975], data.(compDataType).P_AChCBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(compDataType).P_AChCBV_yCI95 = bsxfun(@times, data.(compDataType).P_AChCBV_SEM, data.(compDataType).P_AChCBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(compDataType).P_NECBV_N_Exp = size(data.(compDataType).P_NECBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(compDataType).P_NECBV_Mean = mean(data.(compDataType).P_NECBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(compDataType).P_NECBV_SEM = std(data.(compDataType).P_NECBVRaw,1)/sqrt(data.(compDataType).P_NECBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(compDataType).P_NECBV_CI95 = tinv([0.025 0.975], data.(compDataType).P_NECBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(compDataType).P_NECBV_yCI95 = bsxfun(@times, data.(compDataType).P_NECBV_SEM, data.(compDataType).P_NECBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(compDataType).P_AChGFP_N_Exp = size(data.(compDataType).P_AChGFPRaw,1);
    data.(compDataType).P_AChGFP_Mean = mean(data.(compDataType).P_AChGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(compDataType).P_AChGFP_SEM = std(data.(compDataType).P_AChGFPRaw,1)/sqrt(data.(compDataType).P_AChGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(compDataType).P_AChGFP_CI95 = tinv([0.025 0.975], data.(compDataType).P_AChGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(compDataType).P_AChGFP_yCI95 = bsxfun(@times, data.(compDataType).P_AChGFP_SEM, data.(compDataType).P_AChGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(compDataType).P_NEGFP_N_Exp = size(data.(compDataType).P_NEGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(compDataType).P_NEGFP_Mean = mean(data.(compDataType).P_NEGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(compDataType).P_NEGFP_SEM = std(data.(compDataType).P_NEGFPRaw,1)/sqrt(data.(compDataType).P_NEGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(compDataType).P_NEGFP_CI95 = tinv([0.025 0.975], data.(compDataType).P_NEGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(compDataType).P_NEGFP_yCI95 = bsxfun(@times, data.(compDataType).P_NEGFP_SEM, data.(compDataType).P_NEGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(compDataType).PupilDiameter_N_Exp = size(data.(compDataType).PupilDiameterRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(compDataType).PupilDiameter_Mean = mean(data.(compDataType).PupilDiameterRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(compDataType).PupilDiameter_SEM = std(data.(compDataType).PupilDiameterRaw,1)/sqrt(data.(compDataType).PupilDiameter_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(compDataType).PupilDiameter_CI95 = tinv([0.025 0.975], data.(compDataType).PupilDiameter_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(compDataType).PupilDiameter_yCI95 = bsxfun(@times, data.(compDataType).PupilDiameter_SEM, data.(compDataType).PupilDiameter_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    % data.(compDataType).P_CBV_N_Exp = size(data.(compDataType).P_CBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    % data.(compDataType).P_CBV_Mean = mean(data.(compDataType).P_CBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    % data.(compDataType).P_CBV_SEM = std(data.(compDataType).P_CBVRaw,1)/sqrt(data.(compDataType).P_CBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    % data.(compDataType).P_CBV_CI95 = tinv([0.025 0.975], data.(compDataType).P_CBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    % data.(compDataType).P_CBV_yCI95 = bsxfun(@times, data.(compDataType).P_CBV_SEM, data.(compDataType).P_CBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

end
%% plot the response
summaryFigureN = figure ;
% Subplot 1
% Opto
ax1 = subplot(3,1,1);
p1= plot(data.Opto.mean_timeVector,data.Opto.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on
NEn =  data.Opto.P_AChCBV_Mean + data.Opto.P_AChCBV_yCI95;
patch([data.Opto.mean_timeVector' fliplr(data.Opto.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0.8500 0.3250 0.0980],'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.Opto.mean_timeVector,data.Opto.P_AChCBV_Mean + data.Opto.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)

p2= plot(data.Opto.mean_timeVector,data.Opto.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
NEn =  data.Opto.P_NECBV_Mean + data.Opto.P_NECBV_yCI95;
patch([data.Opto.mean_timeVector' fliplr(data.Opto.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0.6350 0.0780 0.1840],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.Opto.mean_timeVector,data.Opto.P_NECBV_Mean + data.Opto.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)

ylabel('CBV \DeltaF/F (%)')

yyaxis right
p3= plot(data.Opto.mean_timeVector,data.Opto.PupilDiameter_Mean,'-','color','k','LineWidth',2);
hold on
NEn =  data.Opto.PupilDiameter_Mean + data.Opto.PupilDiameter_yCI95;
patch([data.Opto.mean_timeVector' fliplr(data.Opto.mean_timeVector')], [NEn(1,:) fliplr(NEn(2,:))],'k','FaceColor','k','EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.Opto.mean_timeVector,data.Opto.PupilDiameter_Mean + data.Opto.PupilDiameter_yCI95,'-','color','k','LineWidth',0.10)
ylabel('Diameter (Z)')

xline(0,'-',{'Start'});
xline(5,'-',{'End'});
% pV= plot(data.Opto.mean_timeVector,data.Opto.P_CBV_Mean,'-','color','red','LineWidth',2);
% hold on
% plot(data.Opto.mean_timeVector,data.Opto.P_CBV_Mean + data.Contra.P_CBV_yCI95,'-','color','red','LineWidth',0.10)
title('OptoStimulus Evoked Response')
% NE
p4 = plot(data.Opto.mean_timeVector,data.Opto.P_NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(data.Opto.mean_timeVector,data.Opto.P_NEGFP_Mean + data.Opto.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)

%ACh
p5 = plot(data.Opto.mean_timeVector,data.Opto.P_AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
plot(data.Opto.mean_timeVector,data.Opto.P_AChGFP_Mean + data.Opto.P_NEGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
ylabel('\DeltaF/F')
axis square

% ACh - NE
% pq = plot(data.Opto.mean_timeVector, (data.Opto.P_AChGFP_Mean - data.Opto.P_NEGFP_Mean),'-','color',[0 0.2 0.1],'LineWidth',2);

ax1.YAxis(1).Color = 'r';
ax1.YAxis(2).Color = 'k';
xlabel('Peri-stimulus time (s)')
% legend([p1],'CBV LH')
legend([p1,p2,p3,p4,p5],'Contra', 'Ipsi','Pupil','Ipsi GCaMP','Contra GCaMP','Location','best')

% legend([p3,p4,pV],'NE', 'ACh','CBV')
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
% ylim([-0.5 3]);
xlim([-5 15])
%% ECoG Spectrogram LH
ax2 = subplot(3,1,2);
imagesc(data.Opto.mean_T,data.Opto.mean_F,data.Opto.mean_CortS_LH)
title('ECoG LH')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
% yscale log
axis square
axis xy
set(gca,'box','off')
% xlim([-5 15])
xline(0,'-',{'Start'});
xline(5,'-',{'End'});

%% ECoG Spectrogram RH
% ax3 = subplot(3,1,3);
% imagesc(data.Opto.mean_T,data.Opto.mean_F,data.Opto.mean_CortS_RH)
% title('ECoG RH')
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
% linkaxes([ax1,ax2,ax3],'x')
linkaxes([ax1,ax2],'x')
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
% ax3Pos = get(ax3,'position');

ax2Pos(3:4) = ax1Pos(3:4);
% ax3Pos(3:4) = ax1Pos(3:4);
set(ax1,'position',ax1Pos);
set(ax2,'position',ax2Pos);
% set(ax3,'position',ax3Pos);
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
