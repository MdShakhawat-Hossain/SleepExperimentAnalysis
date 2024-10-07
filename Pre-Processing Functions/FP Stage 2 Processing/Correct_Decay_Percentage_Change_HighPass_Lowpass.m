function [expCorrected_percentage,expCorrected_ZScored,Params] = Correct_Decay_Percentage_Change_HighPass_Lowpass(FiberData,Forcesensor,Params, ChannelName,trimStart,trimEnd)

            %% Prefilter the raw to remove high frequency noise
            Params.Data_Freq=12;
            [z,p,k]=butter(3,Params.low_Freq/(0.5*Params.DataFs),'low'); %Low pass for optical data to physiologically relevant range
            [Params.sos_Data,Params.g_Data]=zp2sos(z,p,k);

            RawData=filtfilt(Params.sos_Low,Params.g_Low,FiberData); %Low pass filter data below 1Hz before fitting to remove metabolic clearance/photobleaching trends
            %% remove some data to perform a better fit
            sampleStart = floor(trimStart*Params.DataFs); % samples to be removed from the start
            sampleEnd = floor(trimEnd*Params.DataFs); % samples to be removed from the end
            
            RawData = RawData(sampleStart+1:end-sampleEnd,:); % remove data from the front and end
            Forcesensor = Forcesensor(sampleStart+1:end-sampleEnd); % remove data from the front and end
            figTime=(1:length(RawData(:,1)))/(Params.DataFs*60); % figure time. usefull for ploting the data. 
            %% get the filter properties
            Params.decay_Freq = 0.001;
            Params.high_Freq=0.001;

            [z,p,k]=butter(3,Params.decay_Freq/(0.5*Params.DataFs),'low'); 
            [Params.sos_decay,Params.g_decay]=zp2sos(z,p,k);            

            [z,p,k]=butter(3,Params.high_Freq/(0.5*Params.DataFs),'high'); 
            [Params.sos_high,Params.g_high]=zp2sos(z,p,k);
            %% perform the filters to remove the decay and calculate moving average
            LowPassFitData = filtfilt(Params.sos_decay,Params.g_decay,RawData); % low pass filtered fit line. Also represent the mean. F0. 
            HighPassData = filtfilt(Params.sos_high,Params.g_high,RawData); % High pass filtered data. Also represents the absolute changes. Fluoresence - Mean. Delta F
            %% calculate percentage changes
            expCorrected_percentage.F405 = 100*(HighPassData(:,1)./LowPassFitData(:,1));
            expCorrected_percentage.F465 = 100*(HighPassData(:,2)./LowPassFitData(:,2));
            expCorrected_percentage.F560 = 100*(HighPassData(:,3)./LowPassFitData(:,3));
            %% calculate Z-Scored Changes
            % calculate the standard deviation 
            STD.F405 = sqrt(sum(HighPassData(:,1).^2)/(length(HighPassData(:,1))-1));
            STD.F465 = sqrt(sum(HighPassData(:,2).^2)/(length(HighPassData(:,2))-1));
            STD.F560 = sqrt(sum(HighPassData(:,3).^2)/(length(HighPassData(:,3))-1));

            expCorrected_ZScored.F405 = HighPassData(:,1)./STD.F405;
            expCorrected_ZScored.F465 = HighPassData(:,2)./STD.F465;
            expCorrected_ZScored.F560 = HighPassData(:,3)./STD.F560;
            %% plot the percentage changed data
            %405
            PercentImage405 = figure;

            subplot(811)
            plot(figTime, Forcesensor);            
            ylabel('movement');title('force sensor movement');legend('force sensor')
            set(gca,'Xticklabel',[])

            subplot(8,1,2:3)
            plot(figTime, RawData(:,1));
            hold on
            plot(figTime, LowPassFitData(:,1));
            ylabel('fluoresence value');title('fluoresence decay correction');legend('Raw Data','Fit Line')
            set(gca,'Xticklabel',[])

            subplot(8,1,4:5)
            plot(figTime, HighPassData(:,1));
            set(gca,'Xticklabel',[])
            ylabel('fluoresence value');title('fluoresence decay correction');legend('Corrected Data')

            subplot(8,1,6:8)
            plot(figTime, expCorrected_percentage.F405);
            xlabel('time (min)'); 
            ylabel('\Delta F/F (%)');title('Percentage fluoresence change');

            PercentImage405.WindowState = 'minimized';
            drawnow;
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_percentage' ChannelName '405' '.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_percentage' ChannelName '405' '.tiff'],'tiff')        
            close(PercentImage405)

            % 465
            PercentImage465 = figure;

            subplot(811)
            plot(figTime, Forcesensor);            
            ylabel('movement');title('force sensor movement');legend('force sensor')
            set(gca,'Xticklabel',[])

            subplot(8,1,2:3)
            plot(figTime, RawData(:,2));
            hold on
            plot(figTime, LowPassFitData(:,2));
            ylabel('fluoresence value');title('fluoresence decay correction');legend('Raw Data','Fit Line')
            set(gca,'Xticklabel',[])

            subplot(8,1,4:5)
            plot(figTime, HighPassData(:,2));
            set(gca,'Xticklabel',[])
            ylabel('fluoresence value');title('fluoresence decay correction');legend('Corrected Data')

            subplot(8,1,6:8)
            plot(figTime, expCorrected_percentage.F465);
            xlabel('time (min)'); 
            ylabel('\Delta F/F (%)');title('Percentage fluoresence change');

            PercentImage465.WindowState = 'minimized';
            drawnow;
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_percentage' ChannelName '465' '.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_percentage' ChannelName '465' '.tiff'],'tiff')        
            close(PercentImage465)

            % 560
            PercentImage560 = figure;

            subplot(811)
            plot(figTime, Forcesensor);            
            ylabel('movement');title('force sensor movement');legend('force sensor')
            set(gca,'Xticklabel',[])

            subplot(8,1,2:3)
            plot(figTime, RawData(:,3));
            hold on
            plot(figTime, LowPassFitData(:,3));
            ylabel('fluoresence value');title('fluoresence decay correction');legend('Raw Data','Fit Line')
            set(gca,'Xticklabel',[])

            subplot(8,1,4:5)
            plot(figTime, HighPassData(:,3));
            set(gca,'Xticklabel',[])
            ylabel('fluoresence value');title('fluoresence decay correction');legend('Corrected Data')

            subplot(8,1,6:8)
            plot(figTime, expCorrected_percentage.F560);
            xlabel('time (min)'); 
            ylabel('\Delta F/F (%)');title('Percentage fluoresence change');

            PercentImage560.WindowState = 'minimized';
            drawnow;
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_percentage' ChannelName '560' '.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_percentage' ChannelName '560' '.tiff'],'tiff')        
            close(PercentImage560)
            %% plot the Z Scored data
            % figure;
            % subplot(211)
            % plot(figTime, HighPassData(:,1));
            % hold on
            % plot(figTime, LowPassFitData(:,1));
            % subplot(212)
            % plot(figTime, expCorrected_ZScored.F405);
            % 
            % figure;
            % subplot(211)
            % plot(figTime, HighPassData(:,2));
            % hold on
            % plot(figTime, LowPassFitData(:,2));
            % subplot(212)
            % plot(figTime, expCorrected_ZScored.F465);
            % 
            % figure;
            % subplot(211)
            % plot(figTime, HighPassData(:,3));
            % hold on
            % plot(figTime, LowPassFitData(:,3));
            % subplot(212)
            % plot(figTime, expCorrected_ZScored.F560);
            %% correct motion with 405 relationship
            [expCorrected_percentage] =  control405Correction_Update(expCorrected_percentage,Params,ChannelName);
            [expCorrected_ZScored] =  control405Correction_Update(expCorrected_ZScored,Params,ChannelName);
            %% extract the fit line in 405nm. extract the changes due to motion artifact
            % Params.isos_Freq = 0.002;
            % 
            % [z,p,k]=butter(3,Params.isos_Freq/(0.5*Params.DataFs),'low'); 
            % [Params.sos_isos,Params.g_isos]=zp2sos(z,p,k); 
            % 
            % IsosbesticFitData = filtfilt(Params.sos_isos,Params.g_isos,expCorrected_percentage.F405); % low pass filtered fit line. Also represent the mean. F0. 
            % figure;
            % plot(figTime,expCorrected_percentage.F405); hold on; plot(figTime,IsosbesticFitData);
            % expCorrected_percentage_motion.F405 = expCorrected_percentage.F405 - IsosbesticFitData;
            % plot(figTime,expCorrected_percentage_motion.F405);
            %% correct for large changes in 465 signal due to motion changes in 405
            % expCorrected_percentage_motion.F465 = expCorrected_percentage.F465;% - IsosbesticFitData;
            % figure
            % plot(figTime,expCorrected_percentage_motion.F465);
            % hold on; plot(figTime,IsosbesticFitData);
            % plot(figTime,expCorrected_percentage.F465);
            %% correct for large changes in 560 signal due to motion changes in 465
            % expCorrected_percentage_motion.F560 = expCorrected_percentage.F560; % + IsosbesticFitData;
            % figure
            % plot(figTime,expCorrected_percentage_motion.F560);
            % hold on; plot(figTime,IsosbesticFitData);
            % plot(figTime,expCorrected_percentage.F560);
            %% Remove large spikes in the data. These are usually outliers
            % filt_window = 5;
            % Medfit.F405 = medfilt1(expCorrected_percentage_motion.F405, filt_window*Params.DataFs);
            % Medfit.F465 = medfilt1(expCorrected_percentage_motion.F465, filt_window*Params.DataFs);
            % Medfit.F560 = medfilt1(expCorrected_percentage_motion.F560, filt_window*Params.DataFs);
            % 
            % Residual.F405 = expCorrected_percentage_motion.F405-Medfit.F405;
            % Residual.F465 = expCorrected_percentage_motion.F465-Medfit.F465;
            % Residual.F560 = expCorrected_percentage_motion.F560-Medfit.F560;
            % 
            % spike_thresh = 3;
            % removeSpike.F405 = Residual.F405 >= spike_thresh*std(Medfit.F405);
            % removeSpike.F465 = Residual.F465 >= spike_thresh*std(Medfit.F465);
            % removeSpike.F560 = Residual.F560 >= spike_thresh*std(Medfit.F560);
            % 
            % removeSpike.F405 = Residual.F405 >= spike_thresh*range(Medfit.F405);
            % removeSpike.F465 = Residual.F465 >= spike_thresh*range(Medfit.F465);
            % removeSpike.F560 = Residual.F560 >= spike_thresh*range(Medfit.F560);
            % 
            % if any(removeSpike.F405)
            %     disp('Spikes in 405 channel');
            %     expCorrected_percentage_motion.F405(removeSpike.F405) = Medfit.F405(removeSpike.F405);
            % end
            % 
            % if any(removeSpike.F465)
            %     disp('Spikes in 465 channel');
            %     expCorrected_percentage_motion.F465(removeSpike.F465) = Medfit.F465(removeSpike.F465);
            % end
            % 
            % if any(removeSpike.F560)
            %     disp('Spikes in 560 channel');
            %     expCorrected_percentage_motion.F560(removeSpike.F560) = Medfit.F560(removeSpike.F560);
            % end
            
            % figure;
            % subplot(311)
            % plot(expCorrected_percentage_motion.F405);
            % subplot(312)
            % plot(expCorrected_percentage_motion.F465);
            % subplot(313)
            % plot(expCorrected_percentage_motion.F560);

end

