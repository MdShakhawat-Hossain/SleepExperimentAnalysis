function TDTFiberPhotometry_bilateralGRABNE_Percent(filepath,rawDataFilespath,time_Trim_Start,time_Trim_End)

            if ~isfolder('../Figures/Corrections/')
                mkdir('../Figures/Corrections/')
            end
        %% Check if the raw fiber data exist
        FiberDataFileStruct = dir('*_FiberData.mat');
        FiberDataFiles = {FiberDataFileStruct.name}';

        RAWDataFileStruct = dir('*_RawData.mat');
        RAWDataFiles = {RAWDataFileStruct.name}';
        if size(FiberDataFiles,1) ~= size(RAWDataFiles,1)
            Data = TDTbin2mat(filepath);
            %%
            Params.DataFs=round(Data.streams.x560B.fs);
            Params.RawFs = Data.streams.Fitr.fs;
            Params.low_Freq = 1;
            Params.Fit_Freq=0.05;
            Params.plot_Freq=0.1;
            Params.ball_Freq=10;
            %% Extract the onsets 
            Params.DataOnset = Data.epocs.Top_.onset; % get the onset from the top epoc
            Params.DataOffset = Data.epocs.Bot_.offset;% get the offset from the bot epoc
            %% Filter Parameters
            [z,p,k]=butter(3,Params.low_Freq/(0.5*Params.DataFs),'low'); %Low pass for optical data to physiologically relevant range
            [Params.sos_Low,Params.g_Low]=zp2sos(z,p,k);
        
            [z,p,k]=butter(3,Params.ball_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for locomotion data
            [Params.sos_ball,Params.g_ball]=zp2sos(z,p,k);
        
            [z,p,k]=butter(3,Params.Fit_Freq/(0.5*Params.DataFs),'low'); %design lowpass filter for hemodynamic correction
            [Params.sos_Fit,Params.g_Fit]=zp2sos(z,p,k);
        
            [z,p,k]=butter(3,Params.plot_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for plot
            [Params.sos_plot,Params.g_plot]=zp2sos(z,p,k);            
        end
%%
    for trialNum = 1:1:length(rawDataFilespath)
        %% generate file save path for eACh trials
        % trialNum = 4;
            saveExt = [rawDataFilespath{trialNum}];% 
            strIdx = strfind(saveExt,'RawData');
            Params.savepath =  [saveExt(1:strIdx-1) 'FiberData'];
            fiberFile = [Params.savepath '.mat'];
            %% Fiber signals
            if exist(fiberFile,"file")
                load(fiberFile);
                Raw_LH = FiberData.RawData.Raw_LH;
                Raw_RH = FiberData.RawData.Raw_RH;
                ForceSensor_raw = FiberData.RawData.ForceSensor_raw;
                Params = FiberData.params;
            else
                Raw_LH = [double(Data.streams.x405A.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                    double(Data.streams.x465A.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                    double(Data.streams.x560B.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));]';
                Raw_RH = [double(Data.streams.x405C.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                    double(Data.streams.x465C.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                    double(Data.streams.x560D.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));]';
                
                ForceSensor_raw = double(Data.streams.Forc.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))))';
            end
            %% remove data from the front to fix the exponential decay
            timeN = 15; % remove first 15 seconds;
            sampleStart = floor(timeN*Params.DataFs);
            RawData_LH = Raw_LH;
            RawData_RH = Raw_RH;
            ForceSensor = ForceSensor_raw;
            % adding dummy data to match index
            RawData_LH(1:sampleStart,:) = Raw_LH(sampleStart+1:sampleStart+sampleStart,:);
            RawData_RH(1:sampleStart,:) = Raw_RH(sampleStart+1:sampleStart+sampleStart,:);
            ForceSensor(1:sampleStart,:) = ForceSensor_raw(sampleStart+1:sampleStart+sampleStart,:);

            RawData_LH(sampleStart+1:end,:) = Raw_LH(sampleStart+1:end,:);
            RawData_RH(sampleStart+1:end,:) = Raw_RH(sampleStart+1:end,:);
            ForceSensor(sampleStart+1:end,:) = ForceSensor_raw(sampleStart+1:end,:);
            % add 15s back to first 
            ForceSensor_filt = filtfilt(Params.sos_ball,Params.g_ball,ForceSensor);
            %% flag any weird data from the first hour. Mostly caused by waking up after drug injection        
            if trialNum == 1
                RemovalCheck = figure;
                NPlotTime = (1:length(ForceSensor_filt))./(Params.DataFs); %plot time in minutes
                MPlotTime = (1:length(RawData_LH(:,1)))./(Params.DataFs); %plot time in minutes
                h(1)=subplot(511);
                plot(NPlotTime,ForceSensor_filt)
                xlim([0 NPlotTime(end)]);
                legend('Force Sensor')
                h(2)=subplot(512);
                plot(MPlotTime,RawData_LH(:,2))
                xlim([0 MPlotTime(end)]);
                legend('LH Raw 465')
                title('Raw signals recorded from the fibers')
                h(3)= subplot(513);
                plot(MPlotTime,RawData_LH(:,3))
                xlim([0 MPlotTime(end)]);
                legend('LH Raw 560')
                h(4)=subplot(514);
                plot(MPlotTime,RawData_RH(:,2))
                xlim([0 MPlotTime(end)]);
                legend('RH Raw 465')
                h(5)=subplot(515);
                plot(MPlotTime,RawData_RH(:,3))
                xlim([0 MPlotTime(end)]);
                legend('RH Raw 560')
                ylabel('Raw signals')
                xlabel('Time (min)')
                linkaxes(h,'x');
                drawnow;

                Trimchoice = input('Do you want to remove any part of the data (y/n):  ','s'); disp(' ');
                if Trimchoice == "y"
                    StartDataRemove = input('start time (seconds):  '); disp(' ');
                    EndDataRemove = input('end time (seconds):   '); disp(' ');
                    StartDataSample = (StartDataRemove*Params.DataFs) + 1;
                    EndDataSample = (EndDataRemove*Params.DataFs);
                    % adding dummy data to match index
                    RawDataLH = RawData_LH;
                    RawDataRH = RawData_RH;
                    ForceSensor_F = ForceSensor_filt;

                    RawData_RH = RawDataRH(EndDataSample+1 : end,:);
                    RawData_LH = RawDataLH(EndDataSample+1 : end,:);
                    ForceSensor_filt = ForceSensor_F(EndDataSample+1 : end,:);
                    saveas(gcf,['../Figures/Corrections/' Params.savepath '_WeirdDynamics.fig'],'fig');
                    saveas(gcf,['../Figures/Corrections/' Params.savepath '_WeirdDynamics.tiff'],'tiff'); 
                end
                close(RemovalCheck) 
            end
            %% detect any sudden spikes in the data
            % Outliers are defined as elements more than three 
            % local scaled MAD from the local median over a window length specified by window. This method is also known as a Hampel filter.
            time_N = (1:length(RawData_LH(:,2)))/Params.DataFs;
            [OR_RawData_LH,~] = filloutliers(RawData_LH,"linear","movmedian",round(1*Params.DataFs),"SamplePoints",time_N); % slide using a 0.25 second window
             
            time_N = (1:length(RawData_RH(:,2)))/Params.DataFs;
            [OR_RawData_RH,~] = filloutliers(RawData_RH,"linear","movmedian",round(1*Params.DataFs),"SamplePoints",time_N);
            %% remove outliers using median filter
            [DenoisedData_LH] =  Correct_motionartifact(OR_RawData_LH,Params);  
            [DenoisedData_RH] =  Correct_motionartifact(OR_RawData_RH,Params);
            %% Correct exponential decay
            [expCorrectedData.Percentage.LH,expCorrectedData.ZScored.LH,Params] = Correct_Decay_Percentage_Change_HighPass_Lowpass(DenoisedData_LH,ForceSensor_filt,Params,'ACh',time_Trim_Start,time_Trim_End);
            [expCorrectedData.Percentage.RH,expCorrectedData.ZScored.RH,Params] = Correct_Decay_Percentage_Change_HighPass_Lowpass(DenoisedData_RH,ForceSensor_filt,Params,'NE',time_Trim_Start,time_Trim_End);
            %% patching the removed data
            if trialNum == 1
                if Trimchoice == "y"                    
                    % adding dummy data to percentage and zscored data
                    expCorrectedData_Percentage_LH_F405 = expCorrectedData.Percentage.LH.F405;
                    expCorrectedData_Percentage_RH_F405 = expCorrectedData.Percentage.RH.F405;
                    expCorrectedData_ZScored_LH_F405 = expCorrectedData.ZScored.LH.F405;
                    expCorrectedData_ZScored_RH_F405 = expCorrectedData.ZScored.RH.F405; 

                    expCorrectedData_Percentage_LH_F465 = expCorrectedData.Percentage.LH.F465;
                    expCorrectedData_Percentage_RH_F465 = expCorrectedData.Percentage.RH.F465;
                    expCorrectedData_ZScored_LH_F465 = expCorrectedData.ZScored.LH.F465;
                    expCorrectedData_ZScored_RH_F465 = expCorrectedData.ZScored.RH.F465; 

                    expCorrectedData_Percentage_LH_F560 = expCorrectedData.Percentage.LH.F560;
                    expCorrectedData_Percentage_RH_F560 = expCorrectedData.Percentage.RH.F560;
                    expCorrectedData_ZScored_LH_F560 = expCorrectedData.ZScored.LH.F560;
                    expCorrectedData_ZScored_RH_F560 = expCorrectedData.ZScored.RH.F560; 

                    % add the actual data
                    StremLength = length(expCorrectedData_Percentage_LH_F405);                   
                    
                    expCorrectedData.Percentage.LH.F405(EndDataSample+1:EndDataSample+StremLength) = expCorrectedData_Percentage_LH_F405;
                    expCorrectedData.Percentage.RH.F405(EndDataSample+1:EndDataSample+StremLength)  = expCorrectedData_Percentage_RH_F405;
                    expCorrectedData.ZScored.LH.F405(EndDataSample+1:EndDataSample+StremLength)  = expCorrectedData_ZScored_LH_F405;
                    expCorrectedData.ZScored.RH.F405(EndDataSample+1:EndDataSample+StremLength)  = expCorrectedData_ZScored_RH_F405; 

                    expCorrectedData.Percentage.LH.F465(EndDataSample+1:EndDataSample+StremLength)  = expCorrectedData_Percentage_LH_F465;
                    expCorrectedData.Percentage.RH.F465(EndDataSample+1:EndDataSample+StremLength)  = expCorrectedData_Percentage_RH_F465;
                    expCorrectedData.ZScored.LH.F465(EndDataSample+1:EndDataSample+StremLength)  = expCorrectedData_ZScored_LH_F465;
                    expCorrectedData.ZScored.RH.F465(EndDataSample+1:EndDataSample+StremLength)  = expCorrectedData_ZScored_RH_F465; 

                    expCorrectedData.Percentage.LH.F560(EndDataSample+1:EndDataSample+StremLength)  = expCorrectedData_Percentage_LH_F560;
                    expCorrectedData.Percentage.RH.F560(EndDataSample+1:EndDataSample+StremLength)  = expCorrectedData_Percentage_RH_F560;
                    expCorrectedData.ZScored.LH.F560(EndDataSample+1:EndDataSample+StremLength)  = expCorrectedData_ZScored_LH_F560;
                    expCorrectedData.ZScored.RH.F560(EndDataSample+1:EndDataSample+StremLength)  = expCorrectedData_ZScored_RH_F560; 

                    % add the dummy data in the front
                    expCorrectedData.Percentage.LH.F405(StartDataSample:EndDataSample) = expCorrectedData_Percentage_LH_F405(StartDataSample:EndDataSample);
                    expCorrectedData.Percentage.RH.F405(StartDataSample:EndDataSample)  = expCorrectedData_Percentage_RH_F405(StartDataSample:EndDataSample);
                    expCorrectedData.ZScored.LH.F405(StartDataSample:EndDataSample)  = expCorrectedData_ZScored_LH_F405(StartDataSample:EndDataSample);
                    expCorrectedData.ZScored.RH.F405(StartDataSample:EndDataSample)  = expCorrectedData_ZScored_RH_F405(StartDataSample:EndDataSample); 

                    expCorrectedData.Percentage.LH.F465(StartDataSample:EndDataSample)  = expCorrectedData_Percentage_LH_F465(StartDataSample:EndDataSample);
                    expCorrectedData.Percentage.RH.F465(StartDataSample:EndDataSample)  = expCorrectedData_Percentage_RH_F465(StartDataSample:EndDataSample);
                    expCorrectedData.ZScored.LH.F465(StartDataSample:EndDataSample)  = expCorrectedData_ZScored_LH_F465(StartDataSample:EndDataSample);
                    expCorrectedData.ZScored.RH.F465(StartDataSample:EndDataSample)  = expCorrectedData_ZScored_RH_F465(StartDataSample:EndDataSample); 

                    expCorrectedData.Percentage.LH.F560(StartDataSample:EndDataSample)  = expCorrectedData_Percentage_LH_F560(StartDataSample:EndDataSample);
                    expCorrectedData.Percentage.RH.F560(StartDataSample:EndDataSample)  = expCorrectedData_Percentage_RH_F560(StartDataSample:EndDataSample);
                    expCorrectedData.ZScored.LH.F560(StartDataSample:EndDataSample)  = expCorrectedData_ZScored_LH_F560(StartDataSample:EndDataSample);
                    expCorrectedData.ZScored.RH.F560(StartDataSample:EndDataSample)  = expCorrectedData_ZScored_RH_F560(StartDataSample:EndDataSample);
                    
                    % store the information about weird data removal
                    Params.WData.StartDataRemove = StartDataRemove;
                    Params.WData.EndDataRemove = EndDataRemove;
                    Params.WData.StartDataSample = StartDataSample;
                    Params.WData.EndDataSample = EndDataSample;
                end         
            end            
            %% Converting Tdtomato to CBV
            [AnimalID,~,~] = GetFileInfo_FP(saveExt);
            % some mice were injected with TdTomato instead of a CBV marker. 
            % the actual CBV value is the inverse of tdtomato fluorescence. 
            if AnimalID == "NEACh003"
            expCorrectedData.Percentage.LH.F560 = -expCorrectedData.Percentage.LH.F560;
            expCorrectedData.Percentage.RH.F560 = -expCorrectedData.Percentage.RH.F560;
            expCorrectedData.ZScored.LH.F560 = -expCorrectedData.ZScored.LH.F560;
            expCorrectedData.ZScored.RH.F560 = -expCorrectedData.ZScored.RH.F560; 
            
            elseif AnimalID == "NEACh005"
            expCorrectedData.Percentage.LH.F560 = -expCorrectedData.Percentage.LH.F560;
            expCorrectedData.Percentage.RH.F560 = -expCorrectedData.Percentage.RH.F560;
            expCorrectedData.ZScored.LH.F560 = -expCorrectedData.ZScored.LH.F560;
            expCorrectedData.ZScored.RH.F560 = -expCorrectedData.ZScored.RH.F560;  
            
            elseif AnimalID == "NEACh006"
            expCorrectedData.Percentage.LH.F560 = -expCorrectedData.Percentage.LH.F560;
            expCorrectedData.Percentage.RH.F560 = -expCorrectedData.Percentage.RH.F560;
            expCorrectedData.ZScored.LH.F560 = -expCorrectedData.ZScored.LH.F560;
            expCorrectedData.ZScored.RH.F560 = -expCorrectedData.ZScored.RH.F560;  
            end
            %% hemodynamic Corrections
            CBVCorrectionSavePath = [Params.savepath 'CBVCorrectionParameters' '.mat'];
            if AnimalID == "NEAChM002" || AnimalID ==  "NEAChM003" || AnimalID == "NEACh014" || AnimalID == "NEACh015" || AnimalID == "NEACh016"
                if ~(exist(CBVCorrectionSavePath,'file') == 2)
                [CBVCorrection.Percent.LH.Slope,CBVCorrection.Percent.LH.Inter] = CBVHemodynamicsCorrection(expCorrectedData.Percentage.LH,Params,'Percent','LH');
                [CBVCorrection.Percent.RH.Slope,CBVCorrection.Percent.RH.Inter] = CBVHemodynamicsCorrection(expCorrectedData.Percentage.RH,Params,'Percent','RH');
                [CBVCorrection.ZScored.LH.Slope,CBVCorrection.ZScored.LH.Inter] = CBVHemodynamicsCorrection(expCorrectedData.ZScored.LH,Params,'ZScored','LH');
                [CBVCorrection.ZScored.RH.Slope,CBVCorrection.ZScored.RH.Inter] = CBVHemodynamicsCorrection(expCorrectedData.ZScored.RH,Params,'ZScored','RH');
                save(CBVCorrectionSavePath,"CBVCorrection",'-mat');
                else
                    load(CBVCorrectionSavePath); 
                end
            else % import the parameters used for correction
                if contains(AnimalID,'NE') % these parameters are calculated using GRAB mutants and should only be used for GRAB Sensors
                load('C:\Users\mfh5734\OneDrive - The Pennsylvania State University\Science_Research\SleepExperiment\CBV_Correction_Data\CBVMle_Hemodynamic_Correction_Parameters.mat')
                end
            end
            %% perform the correction
            CBVCorrectionData.Percentage.LH = expCorrectedData.Percentage.LH;
            CBVCorrectionData.Percentage.RH = expCorrectedData.Percentage.RH;

            CBVCorrectionData.Percentage.LH.F465 = CBVCorrectionData.Percentage.LH.F465 - ((CBVCorrectionData.Percentage.LH.F560.*CBVCorrection.Percent.LH.Slope) + CBVCorrection.Percent.LH.Inter);
            CBVCorrectionData.Percentage.RH.F465 = CBVCorrectionData.Percentage.RH.F465 - ((CBVCorrectionData.Percentage.RH.F560.*CBVCorrection.Percent.RH.Slope) + CBVCorrection.Percent.RH.Inter);

            CBVCorrectionData.ZScored.LH = expCorrectedData.ZScored.LH;
            CBVCorrectionData.ZScored.RH = expCorrectedData.ZScored.RH;

            CBVCorrectionData.ZScored.LH.F465 = CBVCorrectionData.ZScored.LH.F465 - ((CBVCorrectionData.ZScored.LH.F560.*CBVCorrection.Percent.LH.Slope) + CBVCorrection.Percent.LH.Inter);
            CBVCorrectionData.ZScored.RH.F465 = CBVCorrectionData.ZScored.RH.F465 - ((CBVCorrectionData.ZScored.RH.F560.*CBVCorrection.Percent.RH.Slope) + CBVCorrection.Percent.RH.Inter);
            %% plot the corrected signal
            plot_CBVCorrectedSignal(expCorrectedData.Percentage,CBVCorrectionData.Percentage,Params,'Percent')    
            plot_CBVCorrectedSignal(expCorrectedData.ZScored,CBVCorrectionData.ZScored,Params,'ZScored')
            %% house keeping
            sampleStart = floor(time_Trim_Start*Params.DataFs);
            sampleEnd= floor(time_Trim_End*Params.DataFs);          
           % Match the data size
            fieldNames = {'F405','F465','F560'};
            for corI=1:1:length(fieldNames)
                    expCorrectedData.Percentage.LH.(fieldNames{corI}) = [expCorrectedData.Percentage.LH.(fieldNames{corI})(1:sampleStart-1); expCorrectedData.Percentage.LH.(fieldNames{corI})(1:end); expCorrectedData.Percentage.LH.(fieldNames{corI})(end-sampleEnd-2:end)];
                    expCorrectedData.Percentage.RH.(fieldNames{corI}) = [expCorrectedData.Percentage.RH.(fieldNames{corI})(1:sampleStart-1); expCorrectedData.Percentage.RH.(fieldNames{corI})(1:end); expCorrectedData.Percentage.RH.(fieldNames{corI})(end-sampleEnd-2:end)];

                    expCorrectedData.ZScored.LH.(fieldNames{corI}) = [expCorrectedData.ZScored.LH.(fieldNames{corI})(1:sampleStart-1); expCorrectedData.ZScored.LH.(fieldNames{corI})(1:end); expCorrectedData.ZScored.LH.(fieldNames{corI})(end-sampleEnd-2:end)];
                    expCorrectedData.ZScored.RH.(fieldNames{corI}) = [expCorrectedData.ZScored.RH.(fieldNames{corI})(1:sampleStart-1); expCorrectedData.ZScored.RH.(fieldNames{corI})(1:end); expCorrectedData.ZScored.RH.(fieldNames{corI})(end-sampleEnd-2:end)];
            end
            
            fieldNames = {'F405','F465','F560'};
            for corI=1:1:length(fieldNames)
                    CBVCorrectionData.ZScored.LH.(fieldNames{corI}) = [CBVCorrectionData.ZScored.LH.(fieldNames{corI})(1:sampleStart-1); CBVCorrectionData.ZScored.LH.(fieldNames{corI})(1:end); CBVCorrectionData.ZScored.LH.(fieldNames{corI})(end-sampleEnd-2:end)];
                    CBVCorrectionData.ZScored.RH.(fieldNames{corI}) = [CBVCorrectionData.ZScored.RH.(fieldNames{corI})(1:sampleStart-1); CBVCorrectionData.ZScored.RH.(fieldNames{corI})(1:end); CBVCorrectionData.ZScored.RH.(fieldNames{corI})(end-sampleEnd-2:end)];
                    
                    CBVCorrectionData.Percentage.LH.(fieldNames{corI}) = [CBVCorrectionData.Percentage.LH.(fieldNames{corI})(1:sampleStart); CBVCorrectionData.Percentage.LH.(fieldNames{corI})(1:end); CBVCorrectionData.Percentage.LH.(fieldNames{corI})(end-sampleEnd-2:end)];
                    CBVCorrectionData.Percentage.RH.(fieldNames{corI}) = [CBVCorrectionData.Percentage.RH.(fieldNames{corI})(1:sampleStart); CBVCorrectionData.Percentage.RH.(fieldNames{corI})(1:end); CBVCorrectionData.Percentage.RH.(fieldNames{corI})(end-sampleEnd-2:end)];
            end
            %% save the data
            FiberData.RawData.Raw_LH = Raw_LH;
            FiberData.RawData.Raw_RH = Raw_RH;
            FiberData.RawData.ForceSensor_raw =  ForceSensor_raw;
                      
            FiberData.pressureSensor  = ForceSensor;
            FiberData.pressureSensor_filt  = ForceSensor_filt;

            FiberData.ACh.expCorrected.ZScored = expCorrectedData.ZScored.LH;
            FiberData.NE.expCorrected.ZScored = expCorrectedData.ZScored.RH;

            FiberData.ACh.dFF0_p = CBVCorrectionData.Percentage.LH;
            FiberData.NE.dFF0_p = CBVCorrectionData.Percentage.RH;

            FiberData.ACh.dFF0_z = CBVCorrectionData.ZScored.LH;
            FiberData.NE.dFF0_z = CBVCorrectionData.ZScored.RH;

            FiberData.params = Params;
            
            clearvars -except -regexp Data$ path$ Params$ trialNum$ time_Trim_Start time_Trim_End
            saveExt = [rawDataFilespath{trialNum}];
            strIdx = strfind(saveExt,'RawData');
            savepathmat =  [saveExt(1:strIdx-1) 'FiberData.mat'];
            save(savepathmat,"FiberData",'-v7.3')
    end