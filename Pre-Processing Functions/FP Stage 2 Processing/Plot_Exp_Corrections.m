% function [expCorrected, predicted,Params] = Correct_ExponentialDecay_HighPass_Update(RawData,forceSensor,Params,ChannelName,timeN)
        % remove some data to perform a better fit
        clear all; close all; clc
        load('NEACh001_230411_13_56_03_FiberData.mat');
            timeN = 210; % remove first 210 seconds;
            Params = FiberData.params;
            sampleN = floor(timeN*Params.DataFs);
            RawData(:,1) = FiberData.NE.rawData.F405;
            RawData(:,2) = FiberData.NE.rawData.F465;
            RawData(:,3) = FiberData.NE.rawData.F560;
            % adding dummy data to match index
            RawData = RawData(sampleN+1:end-sampleN,:);
            Params.decay_Freq=0.001; % 0.001 
            [z,p,k]=butter(3,Params.decay_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for long term decay
            [Params.sos_decay,Params.g_decay]=zp2sos(z,p,k);

            
%% perform exponential correction
        FitData=filtfilt(Params.sos_Fit,Params.g_Fit,RawData); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
        FiltData=filtfilt(Params.sos_Low,Params.g_Low,RawData); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
        
%         speedaccBinary = speedProcess_Doric(forceSensor, Params.DataFs); % for detrend purpose only
%         Params.ExcludeVals = find(speedaccBinary==1); % exclude running data for trend fitting
        %
        % CorrectionChoice = 'y';
        figTime=(1:length(FiltData(:,1)))/(Params.DataFs*60);
%%        
        % while CorrectionChoice == "y"
            % get the filter properties
            [z,p,k]=butter(3,Params.decay_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for long term decay
            [Params.sos_decay,Params.g_decay]=zp2sos(z,p,k);            
            %% Correct Rhodamine blood volume Exp
            DecayImage = figure;
            TESTFItData = filtfilt(Params.sos_decay,Params.g_decay,RawData);
    
            predicted.F560= TESTFItData(:,3);
            expCorrected.F560=FiltData(:,3)-predicted.F560;


            % Plot the exponential fit
            h(1) = subplot(411);
            plot(resample(figTime,10,Params.DataFs),resample(detrend(RawData(:,3),'constant'),10,Params.DataFs)); hold on; plot(resample(figTime,10,Params.DataFs),resample(detrend(predicted.F560,'constant'),10,Params.DataFs)); 
            % title([ 'Detreding of Rhodamine metabolism ' ChannelName]); xlabel('Time (min)'); legend({'Raw Rhodamine brightness','Fit Line'}); 
            xlim([1 figTime(end)]);
            % ylabel('Raw Signal Amplitude');
    
            % h(2)=subplot(412);
            % plot(figTime,detrend(RawData(:,3),'constant')); hold on; plot(figTime,detrend(expCorrected.F560,'constant')); 
            % title(['Detreding of Rhodamine metabolism ' ChannelName]); xlabel('Time (min)'); legend({'Raw Rhodamine brightness','expCorrected Rhodamine'}); xlim([0 figTime(end)]);
            % ylabel('Raw Signal Amplitude');
            %% Correct GRAB-465 Exp
            predicted.F465=TESTFItData(:,2);
            expCorrected.F465=FiltData(:,2)-predicted.F465;
    
            % Plot the exponential fit
            h(3)=subplot(412);
            plot(resample(figTime,10,Params.DataFs),resample(detrend(RawData(:,2),'constant'),10,Params.DataFs)); hold on; plot(resample(figTime,10,Params.DataFs),resample(detrend(predicted.F465,'constant'),10,Params.DataFs)); 

            % plot(figTime,detrend(RawData(:,2),'constant')); hold on; plot(figTime,detrend(predicted.F465,'constant')); plot(figTime,detrend(FitData(:,2),'constant')); 
            % title(['Exponential fit of 465 ' ChannelName]); xlabel('Time (min)'); legend({'Raw 465 brightness','Exponential Fit','Low pass filtered data fit'}); 
            xlim([1 figTime(end)]);
            % ylabel('Raw Signal Amplitude');
            %%
            Params.low_Freq = 0.5;
            [z,p,k]=butter(3,Params.low_Freq/(0.5*Params.DataFs),'low'); %Low pass for optical data to physiologically relevant range
            [Params.sos_Low,Params.g_Low]=zp2sos(z,p,k);
            expCorrected.F465 = filtfilt(Params.sos_Low,Params.g_Low,expCorrected.F465);
            expCorrected.F560 = filtfilt(Params.sos_Low,Params.g_Low,expCorrected.F560);


            h(4)=subplot(413);
            plot(resample(figTime,10,Params.DataFs),resample(detrend(expCorrected.F465,'constant'),10,Params.DataFs),'Color',[0.4660 0.6740 0.1880]); 
            yyaxis("right")
            plot(resample(figTime,10,Params.DataFs),resample(detrend(expCorrected.F560,'constant'),10,Params.DataFs),'Color',[0.6350 0.0780 0.1840]);

            % title(['Exponential fit of 465 ' ChannelName]); xlabel('Time (min)'); 
            legend({'NE','CBV'}); 
            xlim([1 figTime(end)]);
            % ylabel('Raw Signal Amplitude');
            linkaxes(h,'x');
            drawnow;
            %% Correct 405 Exp
            % predicted.F405=TESTFItData(:,1);
            % expCorrected.F405=FiltData(:,1)-predicted.F405;
            % get input from the user about the correction
            % disp(['the highpass frequency is: ' num2str(Params.decay_Freq) 'Hz']);disp(' ')
%             CChoice = input('Is this an accurate exp correction (y/n) : ','s'); disp(' ')
            % CChoice = 'y';
            % if CChoice == "y"
            %     CorrectionChoice = 'n';
            % elseif CChoice == "n"
            %     close(DecayImage);
            %     Params.decay_Freq = input('What is the new decay freq: '); disp(' ')
            % end
        % end
%         DecayImage.WindowState = 'normal';
        % DecayImage.WindowState = 'minimized';
        % drawnow;
        % saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_' ChannelName 'freq_' num2str(Params.decay_Freq)  '.fig'],'fig')
        % saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_' ChannelName 'freq_' num2str(Params.decay_Freq) '.tiff'],'tiff')
        % close(DecayImage)
