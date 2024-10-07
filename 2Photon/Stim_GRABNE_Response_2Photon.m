function [AnalysisResults] = Stim_GRABNE_Response_2Photon(rootFolder,saveFigs,delim,AnalysisResults,FPanimalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S2 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________
% load('S:\NEACh\AnalysisResults_firstHrs.mat');
% AnalysisResults = AnalysisResults_firstHrs;
%% set-up and process data
% FPanimalIDs = {'NEACh005'};
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
        for dd = 1:length(solenoidNames)
            solenoidName = solenoidNames{1,dd};
            data.Z_NE.(solenoidName).count(:,aa) = AnalysisResults.(animalID).Stim.Z_NE.(solenoidName).count;
            % mean
            data.Z_NE.(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).Stim.Z_NE.(solenoidName).GFP.GFP;
            %std
            data.Z_NE.(solenoidName).GFP_std(:,aa) = AnalysisResults.(animalID).Stim.Z_NE.(solenoidName).GFP.GFPStD;
            % neural activity
            % data.cortical.(solenoidName).cortMUA(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).MUA.corticalData;
            % data.cortical.(solenoidName).cortGam(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).Gam.corticalData;
            % data.cortical.(solenoidName).timeVector(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).timeVector;
            % data.cortical.(solenoidName).cortS(:,:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.corticalS;
            % data.cortical.(solenoidName).cortS_Gam(:,:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.corticalS(49:end,20:23);
            % data.cortical.(solenoidName).T(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.T;
            % data.cortical.(solenoidName).F(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.F;
        end
end
% concatenate the data from the contra and ipsi data
% contra
data.Contra.count = data.Z_NE.LPadSol.count;
data.Contra.Z_NEGFP = data.Z_NE.LPadSol.GFP;
data.Contra.Z_NEGFP_std = data.Z_NE.LPadSol.GFP_std;

% data.Contra.cortMUA = data.cortical.LPadSol.cortMUA;
% data.Contra.cortGam = data.cortical.LPadSol.cortGam;
% data.Contra.timeVector = data.cortical.LPadSol.timeVector;
% data.Contra.cortS = data.cortical.LPadSol.cortS;
% data.Contra.cortS_Gam = data.cortical.LPadSol.cortS_Gam;
% data.Contra.T = data.cortical.LPadSol.T;
% data.Contra.F = data.cortical.LPadSol.F;

%ipsi
data.Ipsi.count = data.Z_NE.RPadSol.count;
data.Ipsi.Z_NEGFP_std = data.Z_NE.RPadSol.GFP_std; 
data.Ipsi.Z_NEGFP = data.Z_NE.RPadSol.GFP; 

% data.Ipsi.cortMUA = data.cortical.RPadSol.cortMUA; 
% data.Ipsi.cortGam = data.cortical.RPadSol.cortGam; 
% data.Ipsi.timeVector = data.cortical.RPadSol.timeVector;
% data.Ipsi.cortS = data.cortical.RPadSol.cortS; 
% data.Ipsi.cortS_Gam = data.cortical.RPadSol.cortS_Gam;
% data.Ipsi.T = data.cortical.RPadSol.T;
% data.Ipsi.F = data.cortical.RPadSol.F;

%auditory
data.Auditory.count = data.Z_NE.AudSol.count;

data.Auditory.Z_NEGFP = data.Z_NE.AudSol.GFP;
data.Auditory.Z_NEGFP_std = data.Z_NE.AudSol.GFP_std;

% data.Auditory.cortMUA = data.cortical.AudSol.cortMUA;
% data.Auditory.cortGam = data.cortical.AudSol.cortGam;
% data.Auditory.timeVector = data.cortical.AudSol.timeVector;
% data.Auditory.cortS = data.cortical.AudSol.cortS;
% data.Auditory.cortS_Gam = data.cortical.AudSol.cortS_Gam;
% data.Auditory.T = data.cortical.AudSol.T;
% data.Auditory.F = data.cortical.AudSol.F;


% take the averages of each field through the proper dimension
for ff = 1:length(compDataTypes)
    compDataType = compDataTypes{1,ff};
    data.(compDataType).mean_Count = mean(data.(compDataType).count,2);
    data.(compDataType).std_Count = std(data.(compDataType).count,0,2);   

    data.(compDataType).Z_NEmean_GFP = mean(data.(compDataType).Z_NEGFP,2);
    data.(compDataType).Z_NEstd_GFP = std(data.(compDataType).Z_NEGFP,0,2);

    data.(compDataType).Z_NEmean_GFP_std = mean(data.(compDataType).Z_NEGFP_std,2);

    % 
    % data.(compDataType).mean_CortMUA = mean(data.(compDataType).cortMUA,2);
    % data.(compDataType).std_CortMUA = std(data.(compDataType).cortMUA,0,2);
    % data.(compDataType).mean_CortGam = mean(data.(compDataType).cortGam,2);
    % data.(compDataType).std_CortGam = std(data.(compDataType).cortGam,0,2);
    % data.(compDataType).mean_timeVector = mean(data.(compDataType).timeVector,2);
    % data.(compDataType).mean_CortS = mean(data.(compDataType).cortS,3).*100;
    % data.(compDataType).mean_CortS_Gam = mean(mean(mean(data.(compDataType).cortS_Gam.*100,2),2),3);
    % data.(compDataType).std_CortS_Gam = std(mean(mean(data.(compDataType).cortS_Gam.*100,2),2),0,3);
    % data.(compDataType).mean_T = mean(data.(compDataType).T,2);
    % data.(compDataType).mean_F = mean(data.(compDataType).F,2);
end

%% plot the comparison of GRABNE
summaryFigureN = figure('Name','Stim Evoked');
sgtitle('Stimulus evoked responses in GRABNE signals')

% NE contra stim
ax4 = subplot(1,3,1);
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP + data.Contra.Z_NEmean_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP - data.Contra.Z_NEmean_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
% ax4.YLim = [-6 7];
xlim([-5 15])

% NE ispi stim
ax5 = subplot(1,3,2);
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_GFP + data.Ipsi.Z_NEmean_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_GFP - data.Ipsi.Z_NEmean_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
% ax5.YLim = [-6 7];
xlim([-5 15])

% NE auditory stim
ax6 = subplot(1,3,3);
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_GFP + data.Auditory.Z_NEmean_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_GFP - data.Auditory.Z_NEmean_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
% ax6.YLim = [-6 7];
xlim([-5 15])

% save figure(s)
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim ];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim ];
    end    
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigureN,[dirpath animalID 'Stim-GRABNE']);
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath animalID 'Stim-GRABNE'])
    close 
end

% end
