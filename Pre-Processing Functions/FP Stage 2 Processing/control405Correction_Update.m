function [dF] =  control405Correction_Update(ExpCorrected,Params,ChannelName)
    %% *************
    % NOTE:
    %      In an ideal world, we can use linear fit to perfectly scale the
    %      signals. However, the fit is usually biased by a small portion of
    %      extrem values. Therefore, we need to focus on the majority of the
    %      dataset, instead of all the data, to scale the optical signals

    % use data length to determine how many bins needed
    nbins = floor(length(ExpCorrected.F405)/Params.DataFs/10); % every 10 seconds as one bin
    [N, edges] = histcounts(ExpCorrected.F405, nbins, 'Normalization','probability');
    rareEventIdx = find(N<0.001); % less than 1 percent
    % form an array of idx to remove these rare events
    removeIdx = zeros(length(ExpCorrected.F405),1);
    for idx = 1:length(rareEventIdx)
        edges_idx = rareEventIdx(idx);
        edges_st = edges(edges_idx);
        edges_fin = edges(edges_idx+1);
        remove_tmp = (ExpCorrected.F405>=edges_st)&(ExpCorrected.F405<edges_fin);
        removeIdx = removeIdx+remove_tmp;
    end
    keepIdx = ~removeIdx;
    %% regress FP465 signal against 405 control to get coeffs
    p = polyfit(ExpCorrected.F405(keepIdx), ExpCorrected.F465(keepIdx), 1);
    ScaledData.F465 = polyval(p,ExpCorrected.F405); % this is the signal caused by motion, i.e., independent of NE
    % subtract the scaled signal to get the residual "transients"
    dF.F465 = ExpCorrected.F465- ScaledData.F465; % this is the signal dependent of NE signal

    % ImageCheck = figure;
    % ImageCheck.WindowState = 'minimized';
    % M(1) = subplot(2,1,1);
    % plot(((1:length(dF.F465))/Params.DataFs),dF.F465);
    % hold on 
    % plot(((1:length(ExpCorrected.F465))/Params.DataFs),ExpCorrected.F465);
    % xlim([0 length(ExpCorrected.F465)/Params.DataFs])
    % xlabel('time');
    % ylabel('\deltaF 465');
    % legend('Isosbestic correct','Exp Correct')
    %% regress FP560 signal against 405 control to get coeffs
    p = polyfit(ExpCorrected.F405(keepIdx), ExpCorrected.F560(keepIdx),1);
    ScaledData.F560 = polyval(p,ExpCorrected.F405);%p(1)*ExpCorrected.F405+p(2); % this is the signal caused by motion, i.e., independent of plasma volume
    % subtract the scaled signal to get the residual "transients"
    dF.F560 = ExpCorrected.F560 - ScaledData.F560; % this is the signal dependent of plasma volume signal
    
    % M(2) = subplot(2,1,2);
    % plot(((1:length(dF.F560))/Params.DataFs),dF.F560);
    % hold on 
    % plot(((1:length(ExpCorrected.F560))/Params.DataFs),ExpCorrected.F560);
    % xlim([0 length(ExpCorrected.F560)/Params.DataFs])
    % xlabel('time');
    % ylabel('\deltaF 560');
    % legend('Isosbestic correct','Exp Correct')
    % linkaxes(M);
    %%
    dF.F405 = ExpCorrected.F405;
    % saveas(gcf,['../Figures/Corrections/' Params.savepath 'IsosbesticCorrected_' ChannelName '.fig'],'fig')
    % close(ImageCheck)