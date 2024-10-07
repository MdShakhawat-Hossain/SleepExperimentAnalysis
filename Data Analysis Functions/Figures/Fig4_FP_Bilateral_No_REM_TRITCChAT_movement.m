function [AnalysisResults] = Fig4_FP_Bilateral_No_REM_TRITCChAT_movement(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
IOSanimalIDs = {'T281','T282','T284','T285'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE'};
movecond = {'move_0s','move_1s','move_2s','nomove_0s','nomove_1s','nomove_2s'};
%% FP mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,aa};
    for bb = 1:length(transitions)
        transition = transitions{1,bb};
        for mvn =  1:length(movecond)
        movementname = movecond{1,mvn};
        data.(transition).(movementname).EMG(aa,:) = AnalysisResults.(animalID).Transitions.(transition).(movementname).EMG;
        data.(transition).(movementname).Hip(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).(movementname).Hip;
        data.(transition).(movementname).T = AnalysisResults.(animalID).Transitions.(transition).(movementname).T;
        data.(transition).(movementname).F = AnalysisResults.(animalID).Transitions.(transition).(movementname).F;
        data.(transition).(movementname).Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).(movementname).Cort;
        data.(transition).(movementname).TRITC(aa,:) = AnalysisResults.(animalID).Transitions.(transition).(movementname).TRITC;
        data.(transition).(movementname).GCaMP7s(aa,:) = AnalysisResults.(animalID).Transitions.(transition).(movementname).GCaMP7s;

        data.(transition).(movementname).LH_Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).(movementname).LH_Cort;
        data.(transition).(movementname).LH_TRITC(aa,:) = AnalysisResults.(animalID).Transitions.(transition).(movementname).LH_TRITC;
        data.(transition).(movementname).LH_GCaMP7s(aa,:) = AnalysisResults.(animalID).Transitions.(transition).(movementname).LH_GCaMP7s;
        data.(transition).(movementname).RH_Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).(movementname).RH_Cort;
        data.(transition).(movementname).RH_TRITC(aa,:) = AnalysisResults.(animalID).Transitions.(transition).(movementname).RH_TRITC;
        data.(transition).(movementname).RH_GCaMP7s(aa,:) = AnalysisResults.(animalID).Transitions.(transition).(movementname).RH_GCaMP7s;
        end
    end
end
% take average for each behavioral transition _RH_
for cc = 1:length(transitions)
    transition = transitions{1,cc};
    for mvn =  1:length(movecond)
        movementname = movecond{1,mvn};
        data.(transition).(movementname).meanEMG = nanmean(data.(transition).(movementname).EMG,1);
        data.(transition).(movementname).stdEMG = nanstd(data.(transition).(movementname).EMG,0,1);
        data.(transition).(movementname).meanHip = nanmean(data.(transition).(movementname).Hip,3)*100;
    
        data.(transition).(movementname).meanTRITC = nanmean(data.(transition).(movementname).TRITC,1);
        data.(transition).(movementname).stdTRITC = nanstd(data.(transition).(movementname).TRITC,0,1);
        data.(transition).(movementname).meanGCaMP7s = nanmean(data.(transition).(movementname).GCaMP7s,1);
        data.(transition).(movementname).stdGCaMP7s = nanstd(data.(transition).(movementname).GCaMP7s,0,1);
        data.(transition).(movementname).meanCort = nanmean(data.(transition).(movementname).Cort,3)*100;
    
        data.(transition).(movementname).mean_LH_TRITC = nanmean(data.(transition).(movementname).LH_TRITC,1);
        data.(transition).(movementname).std_LH_TRITC = nanstd(data.(transition).(movementname).LH_TRITC,0,1);
        data.(transition).(movementname).mean_LH_GCaMP7s = nanmean(data.(transition).(movementname).LH_GCaMP7s,1);
        data.(transition).(movementname).std_LH_GCaMP7s = nanstd(data.(transition).(movementname).LH_GCaMP7s,0,1);
        data.(transition).(movementname).mean_LH_Cort = nanmean(data.(transition).(movementname).LH_Cort,3)*100;
    
        data.(transition).(movementname).mean_RH_TRITC = nanmean(data.(transition).(movementname).RH_TRITC,1);
        data.(transition).(movementname).std_RH_TRITC = nanstd(data.(transition).(movementname).RH_TRITC,0,1);
        data.(transition).(movementname).mean_RH_GCaMP7s = nanmean(data.(transition).(movementname).RH_GCaMP7s,1);
        data.(transition).(movementname).std_RH_GCaMP7s = nanstd(data.(transition).(movementname).RH_GCaMP7s,0,1);
        data.(transition).(movementname).mean_RH_Cort = nanmean(data.(transition).(movementname).RH_Cort,3)*100;
    end
end
T1 = -30 + (1/30):(1/30):30;
T2 = -30 + (1/10):(1/10):30;
%% Fig. 4 Plot TRITC Transition with Cortical LFPs
for mvn =  1:length(movecond)
        movementname = movecond{1,mvn};

        summaryFigure_A = figure('Name','Fig4 A');
        sgtitle([movementname '_movement condition'])
        % Awake to NREM
        ax1 = subplot(3,2,1);
        % TRITC and EMG
        p1 = plot(T1,data.AWAKEtoNREM.(movementname).meanTRITC,'-','color',colors('dark candy apple red'),'LineWidth',2);
        hold on
        
        plot(T1,data.AWAKEtoNREM.(movementname).meanTRITC + data.AWAKEtoNREM.(movementname).stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
        plot(T1,data.AWAKEtoNREM.(movementname).meanTRITC - data.AWAKEtoNREM.(movementname).stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
        ylabel('\DeltaRhodamine')
        xlim([-30,30])
        % ylim([-5,50])
        yyaxis right
        
        
        p2 =  plot(T1,data.AWAKEtoNREM.(movementname).meanGCaMP7s,'-','color','k','LineWidth',2);
        hold on
        plot(T1,data.AWAKEtoNREM.(movementname).meanGCaMP7s + data.AWAKEtoNREM.(movementname).stdGCaMP7s,'-','color','k','LineWidth',0.1);
        plot(T1,data.AWAKEtoNREM.(movementname).meanGCaMP7s - data.AWAKEtoNREM.(movementname).stdGCaMP7s,'-','color','k','LineWidth',0.1);
        
        title('Awake to NREM transition')
        xlabel('Time (s)')
        ylabel('\DeltaGCaMP7s','rotation',-90,'VerticalAlignment','bottom')
        
        
        set(gca,'box','off')
        legend([p1,p2],'Rhodamine','GCaMP7s','Location','best','FontSize',5)
        ax1.YAxis(1).Color = colors('dark candy apple red');
        ax1.YAxis(2).Color = 'k';
        % ylim([-1,0.5])
        ax1.TickLength = [0.03,0.03];
        % cort neural
        ax2 = subplot(3,2,3);
        Semilog_ImageSC(T2,data.AWAKEtoNREM.(movementname).F,data.AWAKEtoNREM.(movementname).meanCort,'y')
        axis xy
        c1 = colorbar;
        ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
        caxis([-100,100])
        xlabel('Time (s)')
        ylabel({'Cortical LFP';'Frequency (Hz)'})
        set(gca,'Yticklabel','10^1')
        xlim([-30,30])
        set(gca,'box','off')
        ax2.TickLength = [0.03,0.03];
        % hippocampal neural
        ax3 = subplot(3,2,5);
        Semilog_ImageSC(T2,data.AWAKEtoNREM.(movementname).F,data.AWAKEtoNREM.(movementname).meanHip,'y')
        c2 = colorbar;
        ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
        caxis([-100,100])
        xlabel('Time (s)')
        ylabel({'Hippocampal LFP';'Frequency (Hz)'})
        set(gca,'Yticklabel','10^1')
        xlim([-30,30])
        set(gca,'box','off')
        ax3.TickLength = [0.03,0.03];
        % NREM to Awake
        ax4 = subplot(3,2,2);
        % TRITC and EMG
        plot(T1,data.NREMtoAWAKE.(movementname).meanTRITC,'-','color',colors('dark candy apple red'),'LineWidth',2);
        hold on
        plot(T1,data.NREMtoAWAKE.(movementname).meanTRITC + data.NREMtoAWAKE.(movementname).stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
        plot(T1,data.NREMtoAWAKE.(movementname).meanTRITC - data.NREMtoAWAKE.(movementname).stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
        ylabel('\DeltaRhodamine')
        xlim([-30,30])
        % ylim([-5,50])
        yyaxis right
        
        
        plot(T1,data.NREMtoAWAKE.(movementname).meanGCaMP7s,'-','color','k','LineWidth',2);
        hold on
        plot(T1,data.NREMtoAWAKE.(movementname).meanGCaMP7s + data.NREMtoAWAKE.(movementname).stdGCaMP7s,'-','color','k','LineWidth',0.1);
        plot(T1,data.NREMtoAWAKE.(movementname).meanGCaMP7s - data.NREMtoAWAKE.(movementname).stdGCaMP7s,'-','color','k','LineWidth',0.1);
        
        title('NREM to Awake transition')
        xlabel('Time (s)')
        ylabel('  \DeltaGCaMP7s','rotation',-90,'VerticalAlignment','bottom')
        set(gca,'box','off')
        ax4.YAxis(1).Color = colors('dark candy apple red');
        ax4.YAxis(2).Color = 'k';
        % ylim([-1,0.5])
        ax4.TickLength = [0.03,0.03];
        % cort neural
        ax5 = subplot(3,2,4);
        Semilog_ImageSC(T2,data.NREMtoAWAKE.(movementname).F,data.NREMtoAWAKE.(movementname).meanCort,'y')
        axis xy
        c3 = colorbar;
        ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
        caxis([-100,100])
        xlabel('Time (s)')
        % ylabel({'Cortical LFP';'Frequency (Hz)'})
        set(gca,'Yticklabel','10^1')
        xlim([-30,30])
        set(gca,'box','off')
        ax5.TickLength = [0.03,0.03];
        % hippocampal neural
        ax6 = subplot(3,2,6);
        Semilog_ImageSC(T2,data.NREMtoAWAKE.(movementname).F,data.NREMtoAWAKE.(movementname).meanHip,'y')
        c4 = colorbar;
        ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
        caxis([-100,100])
        xlabel('Time (s)')
        % ylabel({'Hippocampal LFP';'Frequency (Hz)'})
        set(gca,'Yticklabel','10^1')
        xlim([-30,30])
        set(gca,'box','off')
        ax6.TickLength = [0.03,0.03];
        
        %% axes positionns
        ax1Pos = get(ax1,'position');
        ax2Pos = get(ax2,'position');
        ax3Pos = get(ax3,'position');
        ax4Pos = get(ax4,'position');
        ax5Pos = get(ax5,'position');
        ax6Pos = get(ax6,'position');
        
        ax2Pos(3:4) = ax1Pos(3:4);
        ax3Pos(3:4) = ax1Pos(3:4);
        ax4Pos(3:4) = ax1Pos(3:4);
        ax5Pos(3:4) = ax1Pos(3:4);
        ax6Pos(3:4) = ax1Pos(3:4);
        
        set(ax2,'position',ax2Pos);
        set(ax3,'position',ax3Pos);
        set(ax4,'position',ax4Pos);
        set(ax5,'position',ax5Pos);
        set(ax6,'position',ax6Pos);
        %% save figure(s)
        if strcmp(saveFigs,'y') == true
            dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'MovementFree'];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(summaryFigure_A,[dirpath movementname '_movement_condition']);
            set(summaryFigure_A,'PaperPositionMode','auto');
            print('-painters','-dpdf','-fillpage',[dirpath movementname '_movement_condition'])
        end
        close

end
end
