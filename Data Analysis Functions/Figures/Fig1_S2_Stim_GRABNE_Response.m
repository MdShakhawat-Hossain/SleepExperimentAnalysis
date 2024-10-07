function [AnalysisResults] = Fig1_S2_Stim_GRABNE_Response(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FPanimalIDs)
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
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
        for dd = 1:length(solenoidNames)
            solenoidName = solenoidNames{1,dd};
            data.Z_NE.(solenoidName).count(:,aa) = AnalysisResults.(animalID).Stim.Z_NE.(solenoidName).count;
            % mean
            data.Z_NE.(solenoidName).Rhodamine(:,aa) = AnalysisResults.(animalID).Stim.Z_NE.(solenoidName).Rhodamine.Rhodamine;
            data.Z_Ach.(solenoidName).Rhodamine(:,aa) = AnalysisResults.(animalID).Stim.Z_Ach.(solenoidName).Rhodamine.Rhodamine;
            data.Z_NE.(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).Stim.Z_NE.(solenoidName).GFP.GFP;
            data.Z_Ach.(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).Stim.Z_Ach.(solenoidName).GFP.GFP;
            %std
            data.Z_NE.(solenoidName).Rhodamine_std(:,aa) = AnalysisResults.(animalID).Stim.Z_NE.(solenoidName).Rhodamine.RhodamineStD;
            data.Z_Ach.(solenoidName).Rhodamine_std(:,aa) = AnalysisResults.(animalID).Stim.Z_Ach.(solenoidName).Rhodamine.RhodamineStD;
            data.Z_NE.(solenoidName).GFP_std(:,aa) = AnalysisResults.(animalID).Stim.Z_NE.(solenoidName).GFP.GFPStD;
            data.Z_Ach.(solenoidName).GFP_std(:,aa) = AnalysisResults.(animalID).Stim.Z_Ach.(solenoidName).GFP.GFPStD;

            % neural activity
            data.cortical.(solenoidName).cortMUA(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).MUA.corticalData;
            data.cortical.(solenoidName).cortGam(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).Gam.corticalData;
            data.cortical.(solenoidName).timeVector(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).timeVector;
            data.cortical.(solenoidName).cortS(:,:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.corticalS;
            data.cortical.(solenoidName).cortS_Gam(:,:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.corticalS(49:end,20:23);
            data.cortical.(solenoidName).T(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.T;
            data.cortical.(solenoidName).F(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.F;
        end
end
% concatenate the data from the contra and ipsi data
% contra
data.Contra.count = data.Z_NE.LPadSol.count;
data.Contra.Z_AchRhodamine = data.Z_Ach.RPadSol.Rhodamine;
data.Contra.Z_AchGFP = data.Z_Ach.RPadSol.GFP;
data.Contra.Z_NERhodamine = data.Z_NE.LPadSol.Rhodamine;
data.Contra.Z_NEGFP = data.Z_NE.LPadSol.GFP;

data.Contra.Z_AchRhodamine_std = data.Z_Ach.RPadSol.Rhodamine_std;
data.Contra.Z_AchGFP_std = data.Z_Ach.RPadSol.GFP_std;
data.Contra.Z_NERhodamine_std = data.Z_NE.LPadSol.Rhodamine_std;
data.Contra.Z_NEGFP_std = data.Z_NE.LPadSol.GFP_std;

data.Contra.cortMUA = data.cortical.LPadSol.cortMUA;
data.Contra.cortGam = data.cortical.LPadSol.cortGam;
data.Contra.timeVector = data.cortical.LPadSol.timeVector;
data.Contra.cortS = data.cortical.LPadSol.cortS;
data.Contra.cortS_Gam = data.cortical.LPadSol.cortS_Gam;
data.Contra.T = data.cortical.LPadSol.T;
data.Contra.F = data.cortical.LPadSol.F;

%ipsi
data.Ipsi.count = data.Z_NE.RPadSol.count;
data.Ipsi.Z_AchRhodamine_std = data.Z_Ach.LPadSol.Rhodamine_std; 
data.Ipsi.Z_AchGFP_std = data.Z_Ach.LPadSol.GFP_std;
data.Ipsi.Z_NERhodamine_std = data.Z_NE.RPadSol.Rhodamine_std;
data.Ipsi.Z_NEGFP_std = data.Z_NE.RPadSol.GFP_std; 

data.Ipsi.Z_AchRhodamine = data.Z_Ach.LPadSol.Rhodamine; 
data.Ipsi.Z_AchGFP = data.Z_Ach.LPadSol.GFP;
data.Ipsi.Z_NERhodamine = data.Z_NE.RPadSol.Rhodamine;
data.Ipsi.Z_NEGFP = data.Z_NE.RPadSol.GFP; 

data.Ipsi.cortMUA = data.cortical.RPadSol.cortMUA; 
data.Ipsi.cortGam = data.cortical.RPadSol.cortGam; 
data.Ipsi.timeVector = data.cortical.RPadSol.timeVector;
data.Ipsi.cortS = data.cortical.RPadSol.cortS; 
data.Ipsi.cortS_Gam = data.cortical.RPadSol.cortS_Gam;
data.Ipsi.T = data.cortical.RPadSol.T;
data.Ipsi.F = data.cortical.RPadSol.F;

%auditory
data.Auditory.count = data.Z_NE.AudSol.count;

data.Auditory.Z_NERhodamine = data.Z_NE.AudSol.Rhodamine;
data.Auditory.Z_AchRhodamine = data.Z_Ach.AudSol.Rhodamine;
data.Auditory.Z_NEGFP = data.Z_NE.AudSol.GFP;
data.Auditory.Z_AchGFP = data.Z_Ach.AudSol.GFP;

data.Auditory.Z_NERhodamine_std = data.Z_NE.AudSol.Rhodamine_std;
data.Auditory.Z_AchRhodamine_std = data.Z_Ach.AudSol.Rhodamine_std;
data.Auditory.Z_NEGFP_std = data.Z_NE.AudSol.GFP_std;
data.Auditory.Z_AchGFP_std = data.Z_Ach.AudSol.GFP_std;

data.Auditory.cortMUA = data.cortical.AudSol.cortMUA;
data.Auditory.cortGam = data.cortical.AudSol.cortGam;
data.Auditory.timeVector = data.cortical.AudSol.timeVector;
data.Auditory.cortS = data.cortical.AudSol.cortS;
data.Auditory.cortS_Gam = data.cortical.AudSol.cortS_Gam;
data.Auditory.T = data.cortical.AudSol.T;
data.Auditory.F = data.cortical.AudSol.F;

% take the averages of each field through the proper dimension
for ff = 1:length(compDataTypes)
    compDataType = compDataTypes{1,ff};
    data.(compDataType).mean_Count = mean(data.(compDataType).count,2);
    data.(compDataType).std_Count = std(data.(compDataType).count,0,2);   
    
    data.(compDataType).Z_Achmean_Rhodamine = mean(data.(compDataType).Z_AchRhodamine,2);
    data.(compDataType).Z_Achstd_Rhodamine = std(data.(compDataType).Z_AchRhodamine,0,2);

    data.(compDataType).Z_NEmean_Rhodamine = mean(data.(compDataType).Z_NERhodamine,2);
    data.(compDataType).Z_NEstd_Rhodamine = std(data.(compDataType).Z_NERhodamine,0,2);

    data.(compDataType).Z_Achmean_GFP = mean(data.(compDataType).Z_AchGFP,2);
    data.(compDataType).Z_Achstd_GFP = std(data.(compDataType).Z_AchGFP,0,2);

    data.(compDataType).Z_NEmean_GFP = mean(data.(compDataType).Z_NEGFP,2);
    data.(compDataType).Z_NEstd_GFP = std(data.(compDataType).Z_NEGFP,0,2);

    data.(compDataType).Z_Achmean_Rhodamine_std = mean(data.(compDataType).Z_AchRhodamine_std,2);
    data.(compDataType).Z_NEmean_Rhodamine_std = mean(data.(compDataType).Z_NERhodamine_std,2);

    data.(compDataType).Z_Achmean_GFP_std = mean(data.(compDataType).Z_AchGFP_std,2);
    data.(compDataType).Z_NEmean_GFP_std = mean(data.(compDataType).Z_NEGFP_std,2);


    data.(compDataType).mean_CortMUA = mean(data.(compDataType).cortMUA,2);
    data.(compDataType).std_CortMUA = std(data.(compDataType).cortMUA,0,2);
    data.(compDataType).mean_CortGam = mean(data.(compDataType).cortGam,2);
    data.(compDataType).std_CortGam = std(data.(compDataType).cortGam,0,2);
    data.(compDataType).mean_timeVector = mean(data.(compDataType).timeVector,2);
    data.(compDataType).mean_CortS = mean(data.(compDataType).cortS,3).*100;
    data.(compDataType).mean_CortS_Gam = mean(mean(mean(data.(compDataType).cortS_Gam.*100,2),2),3);
    data.(compDataType).std_CortS_Gam = std(mean(mean(data.(compDataType).cortS_Gam.*100,2),2),0,3);
    data.(compDataType).mean_T = mean(data.(compDataType).T,2);
    data.(compDataType).mean_F = mean(data.(compDataType).F,2);
end

%% plot the response
summaryFigureN = figure ;
% Subplot 1
% mSacrlet
ax1 = subplot(2,2,1);
p1= plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_Rhodamine + data.Ipsi.Z_Achstd_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_Rhodamine - data.Ipsi.Z_Achstd_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
title('Stimulus Evoked Response')
ylabel('\DeltaF/F (Z)')
ax1.YLim = [-2 6];

% NE
p2 = plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP + data.Contra.Z_NEstd_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP - data.Contra.Z_NEstd_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)

%ACh
p3 = plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',2);
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_GFP + data.Ipsi.Z_Achstd_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_GFP - data.Ipsi.Z_Achstd_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
ylabel('\DeltaF/F (Z)')
ax1.YAxis(1).Color = 'k';
% ax2.YAxis(2).Color = 'k';
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'CBV','NE', 'ACh')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ax1.YLim = [-2 6];
xlim([-5 15])
axis square

% Subplot 2
% mSacrlet
ax2 = subplot(2,2,2);
p1= plot(data.Ipsi.mean_timeVector,data.Contra.Z_NEmean_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(data.Ipsi.mean_timeVector,data.Contra.Z_NEmean_Rhodamine + data.Contra.Z_NEstd_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Contra.Z_NEmean_Rhodamine - data.Contra.Z_NEstd_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
title('Stimulus Evoked Response')
ylabel('\DeltaF/F (Z)')
ax2.YLim = [-2 6];

% NE
p2 = plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP + data.Contra.Z_NEstd_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP - data.Contra.Z_NEstd_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)

%ACh
p3 = plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',2);
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_GFP + data.Ipsi.Z_Achstd_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_GFP - data.Ipsi.Z_Achstd_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
ylabel('\DeltaF/F (Z)')
ax2.YAxis(1).Color = 'k';
% ax2.YAxis(2).Color = 'k';
xlabel('Peri-stimulus time (s)')
legend([p1,p2,p3],'CBV','NE', 'ACh')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
ax2.YLim = [-2 6];
xlim([-5 15])
axis square
%% save figure(s)
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'lastHrs' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'firstHrs' delim];
    end    
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigureN,[dirpath 'Stim-FiberSignals_response']); %animalID 
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Stim-FiberSignals_response']) % animalID
    close 
end

end
