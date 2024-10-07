function [dFF0_z] =  ZScoreFiberData(Corrected,Params,ChannelName)

    rmStartIndex = Params.DataFs*Params.baselineStartTime;
    rmEndIndex = Params.DataFs*Params.baselineEndTime;
    %% F465
    dF.F465 = Corrected.F465;

    baseline_F465 = dF.F465(rmStartIndex:rmEndIndex);
    dFF0_z.F465 = (dF.F465 - mean(baseline_F465))/std(baseline_F465);
    %% F560
    dF.F560 = Corrected.F560;
   
    baseline_F560 = dF.F560(rmStartIndex:rmEndIndex);
    dFF0_z.F560 = (dF.F560 - mean(baseline_F560))/std(baseline_F560);
    %% F405
    dF.F405 = Corrected.F405;
   
    baseline_F405 = dF.F405(rmStartIndex:rmEndIndex);
    dFF0_z.F405 = (dF.F405 - mean(baseline_F405))/std(baseline_F405);
    %% plot the data
    figTime=(1:length(dFF0_z.F465))/(Params.DataFs*60);
    % figure; 
    h(1) = subplot(311);
    plot(figTime,filtfilt(Params.sos_Low,Params.g_Low,dFF0_z.F560)); 
    title('Exp Corrected Zscored CBV'); xlabel('Time (min)'); xlim([0 figTime(end)]);
    ylabel('ZScored \DeltaF');

    h(2) = subplot(312);
    plot(figTime,filtfilt(Params.sos_Low,Params.g_Low,dFF0_z.F465)); 
    title('Exp Corrected Zscored Green FP'); xlabel('Time (min)');  xlim([0 figTime(end)]);
    ylabel('ZScored \DeltaF');   

    h(3) = subplot(313);
    plot(figTime,filtfilt(Params.sos_Low,Params.g_Low,dFF0_z.F405)); 
    title('Exp Corrected Zscored 405 FP'); xlabel('Time (min)');  xlim([0 figTime(end)]);
    ylabel('ZScored \DeltaF'); 

    linkaxes(h);
    saveas(gcf,['../Figures/Corrections/' Params.savepath 'Zscored_' ChannelName '.fig'],'fig')
    close 