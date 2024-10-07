function plotExpCorrection(fiberData)
        %% plot exponential fit

        % LH
        % Plot the exponential fit for 560
        figure;
        hn(1) = subplot(211);
        plot(figTime,fiberData.rawData.LH.F560); hold on; plot(figTime,predicted.F560);
        title('LH Exponential fit of TRITC metabolism'); xlabel('Time (min)'); legend({'Raw TRITC brightness','Exponential Fit'}); xlim([0 figTime(end)]);
        ylabel('Raw fluoresence change');

        hn(2) = subplot(212);
        plot(figTime,fiberData.rawData.LH.F560); hold on; plot(figTime,expCorrected.F560); 
        title('LH Exponential fit of TRITC metabolism'); xlabel('Time (min)'); legend({'Raw TRITC brightness','expCorrected TRITC'}); xlim([0 figTime(end)]);
        ylabel('Raw fluoresence change');
        linkaxes(hn,'x')

        % Plot the exponential fit for 465
        figure;
        hn(1) = subplot(211);
        plot(figTime,fiberData.rawData.LH.F465); hold on; plot(figTime,predicted.F465); 
        title('LH Exponential fit of 465'); xlabel('Time (min)'); legend({'Raw 465 brightness','Exponential Fit'}); xlim([0 figTime(end)]);
        ylabel('Raw fluoresence change');

        hn(2) = subplot(212);
        plot(figTime,fiberData.rawData.LH.F465); hold on; plot(figTime,expCorrected.F465);
        title('LH Exponential fit of 465'); xlabel('Time (min)'); legend({'Raw 465 brightness','expCorrected 465'}); xlim([0 figTime(end)]);
        ylabel('Raw fluoresence change');
        linkaxes(hn,'x')
        % Plot the exponential fit for 405

        figure; 
        hn(1) = subplot(211);
        plot(figTime,fiberData.rawData.LH.F405); hold on; plot(figTime,predicted.F405); 
        title('LH Exponential fit of 405'); xlabel('Time (min)'); legend({'Raw 405 brightness','Exponential Fit'}); xlim([0 figTime(end)]);
        ylabel('Raw fluoresence change');
    
        hn(2) = subplot(212);
        plot(figTime,fiberData.rawData.LH.F405); hold on; plot(figTime,expCorrected.F405); 
        title('LH Exponential fit of 405'); xlabel('Time (min)'); legend({'Raw 405 brightness','expCorrected 405'}); xlim([0 figTime(end)]);
        ylabel('Raw fluoresence change');
        linkaxes(hn,'x')

        % RH
        % Plot the exponential fit for 560
        figure;
        hn(1) = subplot(211);
        plot(figTime,fiberData.rawData.RH.F560); hold on; plot(figTime,predicted.F560);
        title('Exponential fit of TRITC metabolism'); xlabel('Time (min)'); legend({'Raw TRITC brightness','Exponential Fit'}); xlim([0 figTime(end)]);
        ylabel('Raw fluoresence change');

        hn(2) = subplot(212);
        plot(figTime,fiberData.rawData.RH.F560); hold on; plot(figTime,expCorrected.F560); 
        title('Exponential fit of TRITC metabolism'); xlabel('Time (min)'); legend({'Raw TRITC brightness','expCorrected TRITC'}); xlim([0 figTime(end)]);
        ylabel('Raw fluoresence change');
        linkaxes(hn,'x')

        % Plot the exponential fit for 465
        figure;
        hn(1) = subplot(211);
        plot(figTime,fiberData.rawData.RH.F465); hold on; plot(figTime,predicted.F465); 
        title('Exponential fit of 465'); xlabel('Time (min)'); legend({'Raw 465 brightness','Exponential Fit'}); xlim([0 figTime(end)]);
        ylabel('Raw fluoresence change');

        hn(2) = subplot(212);
        plot(figTime,fiberData.rawData.RH.F465); hold on; plot(figTime,expCorrected.F465);
        title('Exponential fit of 465'); xlabel('Time (min)'); legend({'Raw 465 brightness','expCorrected 465'}); xlim([0 figTime(end)]);
        ylabel('Raw fluoresence change');
        linkaxes(hn,'x')
        % Plot the exponential fit for 405

        figure; 
        hn(1) = subplot(211);
        plot(figTime,fiberData.rawData.RH.F405); hold on; plot(figTime,predicted.F405); 
        title('Exponential fit of 405'); xlabel('Time (min)'); legend({'Raw 405 brightness','Exponential Fit'}); xlim([0 figTime(end)]);
        ylabel('Raw fluoresence change');
    
        hn(2) = subplot(212);
        plot(figTime,fiberData.rawData.RH.F405); hold on; plot(figTime,expCorrected.F405); 
        title('Exponential fit of 405'); xlabel('Time (min)'); legend({'Raw 405 brightness','expCorrected 405'}); xlim([0 figTime(end)]);
        ylabel('Raw fluoresence change');
        linkaxes(hn,'x')
end