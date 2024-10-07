function plot_CBVCorrectedSignal(expCorrectedData,CBVCorrectionData,Params,ScoringType)

            figTime=(1:length(CBVCorrectionData.LH.F465))/(Params.DataFs*60);

            ImageCheck = figure;
            ImageCheck.WindowState = 'minimized';

            k(1) = subplot(2,1,1);
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,expCorrectedData.RH.F465),'b')
            hold on;
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,expCorrectedData.RH.F560),'r')
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,CBVCorrectionData.RH.F465),'c')

            legend('RH 465','RH 560','RH 465C')
            title(['CBV Correction RH ' ScoringType]); xlabel('Time (min)'); xlim([0 figTime(end)]);
            ylabel('\Delta F/F');

            k(2) = subplot(2,1,2);
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,expCorrectedData.LH.F465),'b')
            hold on;
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,expCorrectedData.LH.F560),'r')
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,CBVCorrectionData.LH.F465),'c')

            legend('LH 465','LH 560','LH 465C')
            title(['CBV Correction LH ' ScoringType]); xlabel('Time (min)'); xlim([0 figTime(end)]);
            ylabel('\Delta F/F');

            linkaxes(k,'x');

            saveas(gcf,['../Figures/Corrections/' Params.savepath ScoringType '_CBVCorrected.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath ScoringType '_CBVCorrected.tiff'],'tiff')
            close(ImageCheck)  
end