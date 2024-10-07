function [AnalysisResults] = Fig4_FP_Bilateral_No_REM_RhodamineChAT_movement_compare(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
IOSanimalIDs = {'T282','T285'};%{'T281','T282','T284','T285'};
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
        data.(transition).(movementname).Rhodamine(aa,:) = AnalysisResults.(animalID).Transitions.(transition).(movementname).Rhodamine;
        data.(transition).(movementname).GCaMP7s(aa,:) = AnalysisResults.(animalID).Transitions.(transition).(movementname).GCaMP7s;

        data.(transition).(movementname).LH_Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).(movementname).LH_Cort;
        data.(transition).(movementname).LH_Rhodamine(aa,:) = AnalysisResults.(animalID).Transitions.(transition).(movementname).LH_Rhodamine;
        data.(transition).(movementname).LH_GCaMP7s(aa,:) = AnalysisResults.(animalID).Transitions.(transition).(movementname).LH_GCaMP7s;
        data.(transition).(movementname).RH_Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).(movementname).RH_Cort;
        data.(transition).(movementname).RH_Rhodamine(aa,:) = AnalysisResults.(animalID).Transitions.(transition).(movementname).RH_Rhodamine;
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
    
        data.(transition).(movementname).meanRhodamine = nanmean(data.(transition).(movementname).Rhodamine,1);
        data.(transition).(movementname).stdRhodamine = nanstd(data.(transition).(movementname).Rhodamine,0,1);
        data.(transition).(movementname).meanGCaMP7s = nanmean(data.(transition).(movementname).GCaMP7s,1);
        data.(transition).(movementname).stdGCaMP7s = nanstd(data.(transition).(movementname).GCaMP7s,0,1);
        data.(transition).(movementname).meanCort = nanmean(data.(transition).(movementname).Cort,3)*100;
    
        data.(transition).(movementname).mean_LH_Rhodamine = nanmean(data.(transition).(movementname).LH_Rhodamine,1);
        data.(transition).(movementname).std_LH_Rhodamine = nanstd(data.(transition).(movementname).LH_Rhodamine,0,1);
        data.(transition).(movementname).mean_LH_GCaMP7s = nanmean(data.(transition).(movementname).LH_GCaMP7s,1);
        data.(transition).(movementname).std_LH_GCaMP7s = nanstd(data.(transition).(movementname).LH_GCaMP7s,0,1);
        data.(transition).(movementname).mean_LH_Cort = nanmean(data.(transition).(movementname).LH_Cort,3)*100;
    
        data.(transition).(movementname).mean_RH_Rhodamine = nanmean(data.(transition).(movementname).RH_Rhodamine,1);
        data.(transition).(movementname).std_RH_Rhodamine = nanstd(data.(transition).(movementname).RH_Rhodamine,0,1);
        data.(transition).(movementname).mean_RH_GCaMP7s = nanmean(data.(transition).(movementname).RH_GCaMP7s,1);
        data.(transition).(movementname).std_RH_GCaMP7s = nanstd(data.(transition).(movementname).RH_GCaMP7s,0,1);
        data.(transition).(movementname).mean_RH_Cort = nanmean(data.(transition).(movementname).RH_Cort,3)*100;
    end
end
T1 = -30 + (1/30):(1/30):30;
T2 = -30 + (1/10):(1/10):30;
%% Fig. 4 Plot Rhodamine Transition with Cortical LFPs


        summaryFigure_A = figure('Name','Fig4 A');
        sgtitle([movementname 'movement condition compare'])
        % Awake to NREM
        % 0 movements
        ax1 = subplot(3,2,1);
       
        p1 = plot(T1,data.AWAKEtoNREM.move_0s.meanRhodamine,'-.','color',colors('dark candy apple red'),'LineWidth',1);
        hold on
        
        p2 = plot(T1,data.AWAKEtoNREM.nomove_0s.meanRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',1);
        ylabel('\DeltaRhodamine')
        xlim([-30,30])
        % ylim([-5,50])
        yyaxis right
        
        
        p3 =  plot(T1,data.AWAKEtoNREM.move_0s.meanGCaMP7s,'-.','color','k','LineWidth',1);
        hold on
        p4 = plot(T1,data.AWAKEtoNREM.nomove_0s.meanGCaMP7s,'-','color','k','LineWidth',1);
        
        title('Awake to NREM transition')
        xlabel('Time (s)')
        ylabel('\DeltaGCaMP7s','rotation',-90,'VerticalAlignment','bottom')
        
        
        set(gca,'box','off')
        legend([p1,p2,p3,p4],'move0R','Nmove0R','move0G','Nmove0G','Location','best','FontSize',5)
        ax1.YAxis(1).Color = colors('dark candy apple red');
        ax1.YAxis(2).Color = 'k';
        % ylim([-1,0.5])
        ax1.TickLength = [0.03,0.03];

        % 1s movements
        ax2 = subplot(3,2,3);
       
        p1 = plot(T1,data.AWAKEtoNREM.move_1s.meanRhodamine,'-.','color',colors('dark candy apple red'),'LineWidth',1);
        hold on
        
        p2 = plot(T1,data.AWAKEtoNREM.nomove_1s.meanRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',1);
        ylabel('\DeltaRhodamine')
        xlim([-30,30])
        % ylim([-5,50])
        yyaxis right
        
        
        p3 =  plot(T1,data.AWAKEtoNREM.move_1s.meanGCaMP7s,'-.','color','k','LineWidth',1);
        hold on
        p4 = plot(T1,data.AWAKEtoNREM.nomove_1s.meanGCaMP7s,'-','color','k','LineWidth',1);
        
        title('Awake to NREM transition')
        xlabel('Time (s)')
        ylabel('\DeltaGCaMP7s','rotation',-90,'VerticalAlignment','bottom')
        
        
        set(gca,'box','off')
        legend([p1,p2,p3,p4],'move1R','Nmove1R','move1G','Nmove1G','Location','best','FontSize',5)
        ax2.YAxis(1).Color = colors('dark candy apple red');
        ax2.YAxis(2).Color = 'k';
        % ylim([-1,0.5])
        ax2.TickLength = [0.03,0.03];

         % 2s movements
        ax3 = subplot(3,2,5);
       
        p1 = plot(T1,data.AWAKEtoNREM.move_2s.meanRhodamine,'-.','color',colors('dark candy apple red'),'LineWidth',1);
        hold on
        
        p2 = plot(T1,data.AWAKEtoNREM.nomove_2s.meanRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',1);
        ylabel('\DeltaRhodamine')
        xlim([-30,30])
        % ylim([-5,50])
        yyaxis right

        p3 =  plot(T1,data.AWAKEtoNREM.move_2s.meanGCaMP7s,'-.','color','k','LineWidth',1);
        hold on
        p4 = plot(T1,data.AWAKEtoNREM.nomove_2s.meanGCaMP7s,'-','color','k','LineWidth',1);
        
        title('Awake to NREM transition')
        xlabel('Time (s)')
        ylabel('\DeltaGCaMP7s','rotation',-90,'VerticalAlignment','bottom')
        
        
        set(gca,'box','off')
        legend([p1,p2,p3,p4],'move2R','Nmove2R','move2G','Nmove2G','Location','best','FontSize',5)
        ax3.YAxis(1).Color = colors('dark candy apple red');
        ax3.YAxis(2).Color = 'k';
        % ylim([-1,0.5])
        ax3.TickLength = [0.03,0.03];

        % NREM to Awake

        % 0 movements
        ax4 = subplot(3,2,2);
       
        p1 = plot(T1,data.NREMtoAWAKE.move_0s.meanRhodamine,'-.','color',colors('dark candy apple red'),'LineWidth',1);
        hold on
        
        p2 = plot(T1,data.NREMtoAWAKE.nomove_0s.meanRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',1);
        ylabel('\DeltaRhodamine')
        xlim([-30,30])
        % ylim([-5,50])
        yyaxis right
        
        
        p3 =  plot(T1,data.NREMtoAWAKE.move_0s.meanGCaMP7s,'-.','color','k','LineWidth',1);
        hold on
        p4 = plot(T1,data.NREMtoAWAKE.nomove_0s.meanGCaMP7s,'-','color','k','LineWidth',1);
        
        title('NREM to AWAKE transition')
        xlabel('Time (s)')
        ylabel('\DeltaGCaMP7s','rotation',-90,'VerticalAlignment','bottom')
        
        
        set(gca,'box','off')
        legend([p1,p2,p3,p4],'move0R','Nmove0R','move0G','Nmove0G','Location','best','FontSize',5)
        ax4.YAxis(1).Color = colors('dark candy apple red');
        ax4.YAxis(2).Color = 'k';
        % ylim([-1,0.5])
        ax4.TickLength = [0.03,0.03];

        % 1s movements
        ax5 = subplot(3,2,4);
       
        p1 = plot(T1,data.NREMtoAWAKE.move_1s.meanRhodamine,'-.','color',colors('dark candy apple red'),'LineWidth',1);
        hold on
        
        p2 = plot(T1,data.NREMtoAWAKE.nomove_1s.meanRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',1);
        ylabel('\DeltaRhodamine')
        xlim([-30,30])
        % ylim([-5,50])
        yyaxis right

        p3 =  plot(T1,data.NREMtoAWAKE.move_1s.meanGCaMP7s,'-.','color','k','LineWidth',1);
        hold on
        p4 = plot(T1,data.NREMtoAWAKE.nomove_1s.meanGCaMP7s,'-','color','k','LineWidth',1);
        
        title('NREM to AWAKE transition')
        xlabel('Time (s)')
        ylabel('\DeltaGCaMP7s','rotation',-90,'VerticalAlignment','bottom')

        set(gca,'box','off')
        legend([p1,p2,p3,p4],'move1R','Nmove1R','move1G','Nmove1G','Location','best','FontSize',5)
        ax5.YAxis(1).Color = colors('dark candy apple red');
        ax5.YAxis(2).Color = 'k';
        % ylim([-1,0.5])
        ax5.TickLength = [0.03,0.03];

         % 2s movements
        ax6 = subplot(3,2,6);
       
        p1 = plot(T1,data.NREMtoAWAKE.move_2s.meanRhodamine,'-.','color',colors('dark candy apple red'),'LineWidth',1);
        hold on
        
        p2 = plot(T1,data.NREMtoAWAKE.nomove_2s.meanRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',1);
        ylabel('\DeltaRhodamine')
        xlim([-30,30])
        % ylim([-5,50])
        yyaxis right

        p3 =  plot(T1,data.NREMtoAWAKE.move_2s.meanGCaMP7s,'-.','color','k','LineWidth',1);
        hold on
        p4 = plot(T1,data.NREMtoAWAKE.nomove_2s.meanGCaMP7s,'-','color','k','LineWidth',1);
        
        title('NREM to AWAKE transition')
        xlabel('Time (s)')
        ylabel('\DeltaGCaMP7s','rotation',-90,'VerticalAlignment','bottom')

        set(gca,'box','off')
        legend([p1,p2,p3,p4],'move2R','Nmove2R','move2G','Nmove2G','Location','best','FontSize',5)
        ax6.YAxis(1).Color = colors('dark candy apple red');
        ax6.YAxis(2).Color = 'k';
        % ylim([-1,0.5])
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
            savefig(summaryFigure_A,[dirpath  'movement_condition_compare']);
            set(summaryFigure_A,'PaperPositionMode','auto');
            print('-painters','-dpdf','-fillpage',[dirpath  'movement_condition_compare'])
        end
        close


end
