function [dF,dFF0,dFF0_z] =  control405Correction(CBVCorrected,Params,ChannelName)
    %% *************
    % NOTE:
    %      In an ideal world, we can use linear fit to perfectly scale the
    %      signals. However, the fit is usually biased by a small portion of
    %      extrem values. Therefore, we need to focus on the majority of the
    %      dataset, instead of all the data, to scale the optical signals

    % use data length to determine how many bins needed
    nbins = floor(length(CBVCorrected.F405)/Params.DataFs/10); % every 10 seconds as one bin
    [N, edges] = histcounts(CBVCorrected.F405, nbins, 'Normalization','probability');
    rareEventIdx = find(N<0.001); % less than 1 percent
    % form an array of idx to remove these rare events
    removeIdx = zeros(length(CBVCorrected.F405),1);
    for idx = 1:length(rareEventIdx)
        edges_idx = rareEventIdx(idx);
        edges_st = edges(edges_idx);
        edges_fin = edges(edges_idx+1);
        remove_tmp = (CBVCorrected.F405>=edges_st)&(CBVCorrected.F405<edges_fin);
        removeIdx = removeIdx+remove_tmp;
    end
    keepIdx = ~removeIdx;
    %% regress FP465 signal against 405 control to get coeffs
    p = polyfit(CBVCorrected.F405(keepIdx), CBVCorrected.F465(keepIdx), 1);
    ScaledData.F465 = p(1)*CBVCorrected.F405+p(2); % this is the signal caused by motion, i.e., independent of NE
    % subtract the scaled signal to get the residual "transients"
    dF.F465 = CBVCorrected.F465- ScaledData.F465; % this is the signal dependent of NE signal
    % calculate baseline

    figure;
    M(1) = subplot(2,1,1);
    plot(((1:length(dF.F465))/Params.DataFs),dF.F465);
    hold on 
    plot(((1:length(CBVCorrected.F465))/Params.DataFs),CBVCorrected.F465);
    xlim([0 length(CBVCorrected.F465)/Params.DataFs])
    xlabel('time');
    ylabel('dF 465');
    legend('Final correct','CBV Correct')
    %% write something to calculate the baseline accurately
%     rmTime = 300;
    rmStartIndex = Params.DataFs*Params.baselineStartTime;
    rmEndIndex = Params.DataFs*Params.baselineEndTime;
    baseline_F465 = dF.F465(rmStartIndex:rmEndIndex);
    dFF0.F465 = dF.F465/mean(baseline_F465); % normalize the change to the baseline
%     dFF0_z.F465 = zscore(dF.F465); % z-score
    dFF0_z.F465 = (dF.F465 - mean(baseline_F465))/std(baseline_F465);

%     figure; 
%     subplot(311);plot(dF.F465,'r');
% %     subplot(312);plot(dFF0.F465,'b');
%     subplot(313);plot(dFF0_z.F465,'k');
    %% regress FP565 signal against 405 control to get coeffs
    p = polyfit(CBVCorrected.F405(keepIdx), CBVCorrected.F560(keepIdx), 1);
    ScaledData.F560 = p(1)*CBVCorrected.F405+p(2); % this is the signal caused by motion, i.e., independent of plasma volume
    % subtract the scaled signal to get the residual "transients"
    dF.F560 = CBVCorrected.F560 - ScaledData.F560; % this is the signal dependent of plasma volume signal
    % calculate baseline
    baseline_F560 = dF.F560(rmStartIndex:rmEndIndex);
    dFF0.F560 = dF.F560/mean(baseline_F560);
    dFF0_z.F560 = (dF.F560 - mean(baseline_F560))/std(baseline_F560);

    M(2) = subplot(2,1,2);
    plot(((1:length(dF.F560))/Params.DataFs),dF.F560);
    hold on 
    plot(((1:length(CBVCorrected.F560))/Params.DataFs),CBVCorrected.F560);
    xlim([0 length(CBVCorrected.F560)/Params.DataFs])
    xlabel('time');
    ylabel('dF 465');
    legend('Final correct','CBV Correct')
    linkaxes(M);
    saveas(gcf,['../Figures/Corrections/' Params.savepath 'IsosbesticCorrected_' ChannelName '.fig'],'fig')
    close
    %% plot the data
    figTime=(1:length(dFF0_z.F465))/(Params.DataFs*60);
    figure; 
    h(1) = subplot(211);
    plot(figTime,filtfilt(Params.sos_Low,Params.g_Low,dFF0_z.F560)); 
    title('Isosbestic Corrected Zscored Rhodamine'); xlabel('Time (min)'); xlim([0 figTime(end)]);
    ylabel('ZScored \DeltaF');

    h(2) = subplot(212);
    plot(figTime,filtfilt(Params.sos_Low,Params.g_Low,dFF0_z.F465)); 
    title('Isosbestic Corrected Zscored GFP'); xlabel('Time (min)');  xlim([0 figTime(end)]);
    ylabel('ZScored \DeltaF');   

    linkaxes(h);
    saveas(gcf,['../Figures/Corrections/' Params.savepath 'FinalZscored_' ChannelName '.fig'],'fig')
    close 