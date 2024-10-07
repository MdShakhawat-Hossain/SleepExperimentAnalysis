function [expCorrected, predicted,Params] = Correct_ExponentialDecay_HighPass(RawData,forceSensor,Params,ChannelName)
        % remove some data to perform a better fit
            timeN = 210; % remove first 210 seconds;
            sampleN = floor(timeN*Params.DataFs);
           
            % adding dummy data to match index
            RawData = RawData(sampleN+1:end-sampleN,:);
            forceSensor = forceSensor(sampleN+1:end-sampleN,:);
%% perform exponential correction
        FitData=filtfilt(Params.sos_Fit,Params.g_Fit,RawData); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
        FiltData=filtfilt(Params.sos_Low,Params.g_Low,RawData); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
        
%         speedaccBinary = speedProcess_Doric(forceSensor, Params.DataFs); % for detrend purpose only
%         Params.ExcludeVals = find(speedaccBinary==1); % exclude running data for trend fitting
        %
        CorrectionChoice = 'y';
        figTime=(1:length(FiltData(:,1)))/(Params.DataFs*60);
        
        while CorrectionChoice == "y"
            % get the filter properties
            [z,p,k]=butter(3,Params.decay_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for long term decay
            [Params.sos_decay,Params.g_decay]=zp2sos(z,p,k);            
            %% Correct Rhodamine blood volume Exp
            DecayImage = figure;
            TESTFItData = filtfilt(Params.sos_decay,Params.g_decay,RawData);
    
            predicted.F560= TESTFItData(:,3);
            % Plot the exponential fit
            h(1) = subplot(411);
            plot(figTime,detrend(RawData(:,3),'constant')); hold on; plot(figTime,detrend(predicted.F560,'constant')); plot(figTime,detrend(FitData(:,3),'constant')); %plot(figTime,WheelData);
            title([ 'Detreding of Rhodamine metabolism ' ChannelName]); xlabel('Time (min)'); legend({'Raw Rhodamine brightness','Fit Line','Low pass filtered data fit'}); xlim([0 figTime(end)]);
            ylabel('Raw Signal Amplitude');
    
            expCorrected.F560=FiltData(:,3)-predicted.F560;
    
            h(2)=subplot(412);
            plot(figTime,detrend(RawData(:,3),'constant')); hold on; plot(figTime,detrend(expCorrected.F560,'constant')); 
            title(['Detreding of Rhodamine metabolism ' ChannelName]); xlabel('Time (min)'); legend({'Raw Rhodamine brightness','expCorrected Rhodamine'}); xlim([0 figTime(end)]);
            ylabel('Raw Signal Amplitude');
            %% Correct GRAB-465 Exp
            predicted.F465=TESTFItData(:,2);
    
            % Plot the exponential fit
            h(3)=subplot(413);
            plot(figTime,detrend(RawData(:,2),'constant')); hold on; plot(figTime,detrend(predicted.F465,'constant')); plot(figTime,detrend(FitData(:,2),'constant')); 
            title(['Exponential fit of 465 ' ChannelName]); xlabel('Time (min)'); legend({'Raw 465 brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
            ylabel('Raw Signal Amplitude');
    
            expCorrected.F465=FiltData(:,2)-predicted.F465;
    
            h(4)=subplot(414);
            plot(figTime,detrend(RawData(:,2),'constant')); hold on; plot(figTime,detrend(expCorrected.F465,'constant')); 
            title(['Exponential fit of 465 ' ChannelName]); xlabel('Time (min)'); legend({'Raw 465 brightness','expCorrected 465'}); xlim([0 figTime(end)]);
            ylabel('Raw Signal Amplitude');
            linkaxes(h,'x');
            DecayImage.WindowState = 'maximized';
            drawnow;
            %% Correct 405 Exp
            predicted.F405=TESTFItData(:,1);
            expCorrected.F405=FiltData(:,1)-predicted.F405;
            % get input from the user about the correction
            disp(['the highpass frequency is: ' num2str(Params.decay_Freq) 'Hz']);disp(' ')
            CChoice = input('Is this an accurate exp correction (y/n) : ','s'); disp(' ')
            if CChoice == "y"
                CorrectionChoice = 'n';
            elseif CChoice == "n"
                close(DecayImage);
                Params.decay_Freq = input('What is the new decay freq: '); disp(' ')
            end
        end
        DecayImage.WindowState = 'normal';
        DecayImage.WindowState = 'minimized';
        drawnow;
        saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_' ChannelName 'freq_' num2str(Params.decay_Freq)  '.fig'],'fig')
        saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_' ChannelName 'freq_' num2str(Params.decay_Freq) '.tiff'],'tiff')
        close(DecayImage)
