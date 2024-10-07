function [expCorrected, predicted,Params] = Correct_ExponentialDecay_HighPass_Update(RawData,forceSensor,Params,ChannelName,timeN,saveChoice)
        % remove some data to perform a better fit
            % timeN = 210; % remove first 210 seconds;
            sampleN = floor(timeN*Params.DataFs);
           
            % adding dummy data to match index
            RawData = RawData(sampleN+1:end-sampleN,:);
            forceSensor = forceSensor(sampleN+1:end-sampleN,:);
%% perform exponential correction
        % FitData=filtfilt(Params.sos_Fit,Params.g_Fit,RawData); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
        FiltData=filtfilt(Params.sos_Low,Params.g_Low,RawData); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
        
%         speedaccBinary = speedProcess_Doric(forceSensor, Params.DataFs); % for detrend purpose only
%         Params.ExcludeVals = find(speedaccBinary==1); % exclude running data for trend fitting
        %
        CorrectionChoice = 'y';
        figTime=(1:length(RawData(:,1)))/(Params.DataFs*60);
        
        while CorrectionChoice == "y"
            % get the filter properties
            [z,p,k]=butter(3,Params.decay_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for long term decay
            [Params.sos_decay,Params.g_decay]=zp2sos(z,p,k);            
            %% Correct Rhodamine blood volume Exp
            DecayImage = figure;
            TESTFitData = filtfilt(Params.sos_decay,Params.g_decay,RawData);

            %%
            Params.high_Freq=0.001; % 0.001

            [z,p,k]=butter(3,Params.high_Freq/(0.5*Params.DataFs),'high'); %Low pass filter for long term decay
            [Params.sos_high,Params.g_high]=zp2sos(z,p,k);

            TESTHighData = filtfilt(Params.sos_high,Params.g_high,RawData);

            figure;
            subplot(311)
            plot(figTime, TESTHighData(:,3));
            subplot(312)
            plot(figTime, TESTFitData(:,3));
            subplot(313)
            plot(figTime, 100*(TESTHighData(:,3)./TESTFitData(:,3)));

            figure;
            subplot(311)
            plot(figTime, TESTHighData(:,2));
            subplot(312)
            plot(figTime, TESTFitData(:,2));
            subplot(313)
            plot(figTime, 100*(TESTHighData(:,2)./TESTFitData(:,2)));

            %%
            predicted.F560= TESTFitData(:,3);
            % Plot the exponential fit
            h(1) = subplot(611);
            
            plot(figTime,detrend(RawData(:,3),'constant')); hold on; plot(figTime,detrend(predicted.F560,'constant')); plot(figTime,detrend(FitData(:,3),'constant')); %plot(figTime,WheelData);
            title([ 'Detreding of CBV metabolism ' ChannelName]); xlabel('Time (min)'); legend({'Raw CBV brightness','Fit Line','Low pass filtered data fit'}); xlim([0 figTime(end)]);
            ylabel('Raw Signal Amplitude');
    
            expCorrected.F560=FiltData(:,3)-predicted.F560;
    
            h(2)=subplot(612);
            plot(figTime,detrend(RawData(:,3),'constant')); hold on; plot(figTime,detrend(expCorrected.F560,'constant')); 
            title(['Detreding of CBV metabolism ' ChannelName]); xlabel('Time (min)'); legend({'Raw CBV brightness','expCorrected Rhodamine'}); xlim([0 figTime(end)]);
            ylabel('Raw Signal Amplitude');
            %% Correct GRAB-465 Exp
            predicted.F465=TESTFitData(:,2);
    
            % Plot the exponential fit
            h(3)=subplot(613);
            plot(figTime,detrend(RawData(:,2),'constant')); hold on; plot(figTime,detrend(predicted.F465,'constant')); plot(figTime,detrend(FitData(:,2),'constant')); 
            title(['Exponential fit of 465 ' ChannelName]); xlabel('Time (min)'); legend({'Raw 465 brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
            ylabel('Raw Signal Amplitude');
    
            expCorrected.F465=FiltData(:,2)-predicted.F465;
    
            h(4)=subplot(614);
            plot(figTime,detrend(RawData(:,2),'constant')); hold on; plot(figTime,detrend(expCorrected.F465,'constant')); 
            title(['Exponential fit of 465 ' ChannelName]); xlabel('Time (min)'); legend({'Raw 465 brightness','expCorrected 465'}); xlim([0 figTime(end)]);
            ylabel('Raw Signal Amplitude');
%             DecayImage.WindowState = 'maximized';
           
            %% Correct 405 Exp
            predicted.F405=TESTFitData(:,1);
            % Plot the exponential fit
            h(5)=subplot(615);
            plot(figTime,detrend(RawData(:,1),'constant')); hold on; plot(figTime,detrend(predicted.F405,'constant')); plot(figTime,detrend(FitData(:,2),'constant')); 
            title(['Exponential fit of 405 ' ChannelName]); xlabel('Time (min)'); legend({'Raw 405 brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
            ylabel('Raw Signal Amplitude');
    
            expCorrected.F405=FiltData(:,1)-predicted.F405;
    
            h(6)=subplot(616);
            plot(figTime,detrend(RawData(:,1),'constant')); hold on; plot(figTime,detrend(expCorrected.F405,'constant')); 
            title(['Exponential fit of 405 ' ChannelName]); xlabel('Time (min)'); legend({'Raw 405 brightness','expCorrected 405'}); xlim([0 figTime(end)]);
            ylabel('Raw Signal Amplitude');
            linkaxes(h,'x');
            drawnow;
            %%
            % get input from the user about the correction
            disp(['the highpass frequency is: ' num2str(Params.decay_Freq) 'Hz']);disp(' ')
%             CChoice = input('Is this an accurate exp correction (y/n) : ','s'); disp(' ')
            CChoice = 'y';
            if CChoice == "y"
                CorrectionChoice = 'n';
            elseif CChoice == "n"
                close(DecayImage);
                Params.decay_Freq = input('What is the new decay freq: '); disp(' ')
            end
        end
%         DecayImage.WindowState = 'normal';
        DecayImage.WindowState = 'minimized';
        drawnow;
        if saveChoice == 'y'
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_' ChannelName 'freq_' num2str(Params.decay_Freq)  '.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_' ChannelName 'freq_' num2str(Params.decay_Freq) '.tiff'],'tiff')
        elseif saveChoice == 'n'
         
        end
        close(DecayImage)
        close all
        %% final expcorrected signal
        expCorrected.F405 = detrend(expCorrected.F405,'constant');
        expCorrected.F465 = detrend(expCorrected.F465,'constant');
        expCorrected.F560 = detrend(expCorrected.F560,'constant');
        %% plot the signal
            % figure;
            subplot(211)
            plot(figTime,expCorrected.F465); ylabel('Percentage Fluoresence Change');
            yyaxis right

            plot(figTime,expCorrected.F560); ylabel('Percentage Fluoresence Change');
            title(['Exp corrected 465vs560 ' ChannelName]); xlabel('Time (min)'); legend({'465','560'}); xlim([0 figTime(end)]);
            drawnow;
            if saveChoice == 'y'
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_' ChannelName '465vs560.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_' ChannelName '465vs560.tiff'],'tiff')
         
            end
            close all

