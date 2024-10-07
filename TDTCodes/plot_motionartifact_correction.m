    function plot_motionartifact_correction(fiberData)

    figure;
    figtime = (1:length(expCorrected.F405))./(Params.DataFs*60);
    h(1) = subplot(311);
    plot(figtime,expCorrected.F405,'k'); hold on; plot(figtime,Medfit.F405,'r');
    title('Motion artifact removed 405'); ylabel('\DeltaF/F');
    legend('expCorrected 405','Median Filtered 405');
    h(2) = subplot(312);
    plot(figtime,expCorrected.F465,'k'); hold on; plot(figtime,Medfit.F465,'r');
    title('Motion artifact removed 465');ylabel('\DeltaF/F');
    legend('expCorrected 465','Median Filtered 465');
    h(3) = subplot(313);
    plot(figtime,expCorrected.F560,'k'); hold on; plot(figtime,Medfit.F560,'r');
    title('Motion artifact removed 560');ylabel('\DeltaF/F');xlabel('Time(min)')
    legend('expCorrected 560','Median Filtered 560');
    linkaxes(h,'x');