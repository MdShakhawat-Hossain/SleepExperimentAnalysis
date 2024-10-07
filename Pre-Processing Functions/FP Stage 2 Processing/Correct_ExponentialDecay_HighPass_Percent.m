function [expCorrected, predicted,Params] = Correct_ExponentialDecay_HighPass_Percent(RawData,forceSensor,Params,ChannelName,timeN,saveChoice)
        % remove some data to perform a better fit
            % timeN = 210; % remove first 210 seconds;
            sampleN = floor(timeN*Params.DataFs);
           
            % adding dummy data to match index
            RawData = RawData(sampleN+1:end-sampleN,:);
            % forceSensor = forceSensor(sampleN+1:end-sampleN,:);
%% perform exponential correction
        % FitData=filtfilt(Params.sos_Fit,Params.g_Fit,RawData); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
        % FiltData=filtfilt(Params.sos_Low,Params.g_Low,RawData); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
        
%         speedaccBinary = speedProcess_Doric(forceSensor, Params.DataFs); % for detrend purpose only
%         Params.ExcludeVals = find(speedaccBinary==1); % exclude running data for trend fitting
        %
        % figTime=(1:length(FiltData(:,1)))/(Params.DataFs*60);
        
            % get the filter properties
            % [z,p,k]=butter(3,Params.decay_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for long term decay
            % [Params.sos_decay,Params.g_decay]=zp2sos(z,p,k);            
            %% Correct Rhodamine blood volume Exp
            % DecayImage = figure;
            % TESTFItData = filtfilt(Params.sos_decay,Params.g_decay,RawData);

            %% calculate the percentage
                DataLength = length(RawData(:,1));
                DataSize =  Params.DataFs*60*2; % 5 minutes bin
                DataHolder = 0;
                Round_Check = 0;
            while 1
                % if Dn == 1
                %     DataHolder = 0;
                % end
                % DataStart = 1 + DataHolder ;
                % DataEnd = DataStart + DataSize;
                % if DataEnd > DataLength
                %     DataEnd = DataLength;
                % end
                % DataHolder = DataEnd;

                 if DataHolder == 0
                        
                        DataStart = 1 + DataHolder;
                        DataEnd = DataStart + DataSize;
                 elseif DataHolder > 0
                        DataStart = DataHolder + 100;%round(15*DataSize/16);
                        DataEnd = DataStart + DataSize;
                 end

                if DataEnd > DataLength
                    DataEnd = DataLength;
                end

                % if Round_Check == 0
                    % DataHolder = DataEnd;
                % else
                    DataHolder = DataStart;
                % end




                % F405
                dFF0_percentage.F405(DataStart:DataEnd) =  (((RawData(DataStart:DataEnd,1)-mean(RawData(DataStart:DataEnd,1))))./mean(RawData(DataStart:DataEnd,1)));
                % F465
                % temp_F465 = dF.F465(DataStart:DataEnd)./dF.F405(DataStart:DataEnd);
        
                % dFF0_percentage.F465(DataStart:DataEnd) = (100.*(temp_F465-mean(temp_F465)))./mean(temp_F465);
                dFF0_percentage.F465(DataStart:DataEnd) = (((RawData(DataStart:DataEnd,2)-mean(RawData(DataStart:DataEnd,2))))./mean(RawData(DataStart:DataEnd,2)));
        
                % F560    
                dFF0_percentage.F560(DataStart:DataEnd) = (((RawData(DataStart:DataEnd,3)-mean(RawData(DataStart:DataEnd,3))))./mean(RawData(DataStart:DataEnd,3)));
        
                    if DataEnd >= DataLength
                        break
                    end

                Round_Check = Round_Check+1;
            end
    %% plot the data
    figTime=(1:length(dFF0_percentage.F465))/(Params.DataFs*60);
    figure; 
    h(1) = subplot(311);
    plot(figTime,filtfilt(Params.sos_Low,Params.g_Low,dFF0_percentage.F560)); 
    title('Exp Corrected percentagescored CBV'); xlabel('Time (min)'); xlim([0 figTime(end)]);
    ylabel('percentageScored \DeltaF');

    h(2) = subplot(312);
    plot(figTime,filtfilt(Params.sos_Low,Params.g_Low,dFF0_percentage.F465)); 
    title('Exp Corrected percentagescored Green FP'); xlabel('Time (min)');  xlim([0 figTime(end)]);
    ylabel('percentageScored \DeltaF');   

    h(3) = subplot(313);
    plot(figTime,filtfilt(Params.sos_Low,Params.g_Low,dFF0_percentage.F405)); 
    title('Exp Corrected percentagescored 405 FP'); xlabel('Time (min)');  xlim([0 figTime(end)]);
    ylabel('percentageScored \DeltaF'); 
%%
    linkaxes(h);
    saveas(gcf,['../Figures/Corrections/' Params.savepath 'Percentage_' ChannelName '.fig'],'fig')
    close
        
%         DecayImage.WindowState = 'normal';
        % DecayImage.WindowState = 'minimized';
        % drawnow;
        % if saveChoice == 'y'
        %     saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_' ChannelName 'freq_' num2str(Params.decay_Freq)  '.fig'],'fig')
        %     saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_' ChannelName 'freq_' num2str(Params.decay_Freq) '.tiff'],'tiff')
        % elseif saveChoice == 'n'
        % 
        % end
        % close(DecayImage)
