
    function plot_405corrections()
%% plot the figures
    figTime=(1:length(RawData(:,1)))/(Params.DataFs*60);
    
    figure;
    h(1) = subplot(611);
    plot(figTime,AnalogWheel);
    title('Analog Wheel/Force Sensor'); ylabel('Activity')
    h(2) = subplot(612);
    plot(figTime, RawData(:,1), 'k'); hold on; plot(figTime, RawData(:,2), 'b'); plot(figTime, RawData(:,3), 'r');
    legend({'405: isobestic','465', '560'});
    ylabel('Raw Signal');
    title('Fiber Photometry Data')
    h(3) = subplot(613);
    plot(figTime, expCorrected.F405, 'k'); hold on; plot(figTime, expCorrected.F465, 'b'); plot(figTime, expCorrected.F560, 'r');
    legend({'405: isobestic','465', '560'});
    ylabel('Detrended Signal');
    h(4) = subplot(614);
    plot(figTime,dF.F465,'k'); hold on;
    plot(figTime,smooth(dF.F465,0.005,'lowess'),'r');
    legend('465', 'Lowpass filtered')
    ylabel('\DeltaF/F');
    h(5) = subplot(615);
    plot(figTime,dF.F560,'k'); hold on;
    plot(figTime,smooth(dF.F560,0.005,'lowess'),'r');
    legend('CBV', 'Lowpass filtered')
    ylabel('\DeltaF/F');
    xlabel('Time (min)');
    
    linkaxes(h,'x');