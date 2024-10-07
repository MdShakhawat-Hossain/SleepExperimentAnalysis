function [OptoStimData] = OptoStimDataOrganize(NECBVData, NEGFPData, AChCBVData, AChStimGFPData, StimPupilData, StimCortMUAData_LH, StimCortGamData_LH, StimCortMUAData_RH, StimCortGamData_RH, StimFileIDs, StimFileEventTimes,solenoid,trialDuration_sec,offset,samplingRate,animalID,AllSpecData, specSamplingRate,timeVector)
        ii = 1;
                    for hh = 1:size(NECBVData,1)
                        stimStartTime = round(StimFileEventTimes(hh,1),1) - 3;
                        stimEndTime = stimStartTime + 15;
                        finalOptoStimFileID = StimFileIDs{hh,1};
                        if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                            AChstimCBVarray = AChCBVData(hh,:);
                            NEstimCBVarray = NECBVData(hh,:);
                            AChstimGFParray = AChStimGFPData(hh,:);
                            NEstimGFParray = NEGFPData(hh,:);

                            stimCortMUAarray_LH = StimCortMUAData_LH(hh,:);
                            stimCortGamArray_LH = StimCortGamData_LH(hh,:);
                            stimCortMUAarray_RH = StimCortMUAData_RH(hh,:);
                            stimCortGamArray_RH = StimCortGamData_RH(hh,:);
        
                            stimPupilarray = StimPupilData(hh,:);
            
                            AChfiltOptoStimCBVarray = sgolayfilt(AChstimCBVarray,3,9);
                            NEfiltOptoStimCBVarray = sgolayfilt(NEstimCBVarray,3,9);
                            AChfiltOptoStimGFParray = sgolayfilt(AChstimGFParray,3,9);
                            NEfiltOptoStimGFParray = sgolayfilt(NEstimGFParray,3,9);
                            PupilfiltOptoStimarray = sgolayfilt(stimPupilarray,3,9);
            
                            AChprocOptoStimCBVData(hh,:) = AChfiltOptoStimCBVarray - mean(AChfiltOptoStimCBVarray(1:(offset*samplingRate)));
                            NEprocOptoStimCBVData(hh,:) = NEfiltOptoStimCBVarray - mean(NEfiltOptoStimCBVarray(1:(offset*samplingRate)));
                            AChprocOptoStimGFPData(hh,:) = AChfiltOptoStimGFParray - mean(AChfiltOptoStimGFParray(1:(offset*samplingRate)));
                            NEprocOptoStimGFPData(hh,:) = NEfiltOptoStimGFParray - mean(NEfiltOptoStimGFParray(1:(offset*samplingRate)));
                            procOptoStimPupilData(hh,:) = PupilfiltOptoStimarray - mean(PupilfiltOptoStimarray(1:(offset*samplingRate)));
          
                            procOptoStimCortMUAData_LH(hh,:) = stimCortMUAarray_LH - mean(stimCortMUAarray_LH(1:(offset*samplingRate)));
                            procOptoStimCortGamData_LH(hh,:) = stimCortGamArray_LH - mean(stimCortGamArray_LH(1:(offset*samplingRate)));
                            
                            procOptoStimCortMUAData_RH(hh,:) = stimCortMUAarray_RH - mean(stimCortMUAarray_RH(1:(offset*samplingRate)));
                            procOptoStimCortGamData_RH(hh,:) = stimCortGamArray_RH - mean(stimCortGamArray_RH(1:(offset*samplingRate)));
                            
                            finalOptoStimStartTimes(ii,1) = stimStartTime;
                            finalOptoStimEndTimes(ii,1) = stimEndTime;
                            finalOptoStimFiles{ii,1} = finalOptoStimFileID;
                            ii = ii + 1;
                        end
                    end

                    %% avearge and standard deviation
                    AChmeanOptoStimCBVData = mean(AChprocOptoStimCBVData,1);
                    AChstdOptoStimCBVData = std(AChprocOptoStimCBVData,0,1);
                    NEmeanOptoStimCBVData = mean(NEprocOptoStimCBVData,1);
                    NEstdOptoStimCBVData = std(NEprocOptoStimCBVData,0,1);
                    AChmeanOptoStimGFPData = mean(AChprocOptoStimGFPData,1);
                    AChstdOptoStimGFPData = std(AChprocOptoStimGFPData,0,1);
                    NEmeanOptoStimGFPData = mean(NEprocOptoStimGFPData,1);
                    NEstdOptoStimGFPData = std(NEprocOptoStimGFPData,0,1);
            
                    meanOptoStimCortMUAData_LH = mean(procOptoStimCortMUAData_LH,1)*100;
                    stdOptoStimCortMUAData_LH = std(procOptoStimCortMUAData_LH,0,1)*100;
                    meanOptoStimCortGamData_LH = mean(procOptoStimCortGamData_LH,1)*100;
                    stdOptoStimCortGamData_LH = std(procOptoStimCortGamData_LH,0,1)*100;
        
                    meanOptoStimCortMUAData_RH = mean(procOptoStimCortMUAData_RH,1)*100;
                    stdOptoStimCortMUAData_RH = std(procOptoStimCortMUAData_RH,0,1)*100;
                    meanOptoStimCortGamData_RH = mean(procOptoStimCortGamData_RH,1)*100;
                    stdOptoStimCortGamData_RH = std(procOptoStimCortGamData_RH,0,1)*100;
        
                    meanOptoStimPupilData = mean(procOptoStimPupilData,1);
                    stdOptoStimPupilData = std(procOptoStimPupilData,0,1);
                    %% extract ECoG spectrograms associated with the stimuli indecies
                    % left hemispheres
                    stimCortZhold_LH = [];
                    for jj = 1:length(finalOptoStimFiles)
                        % load normalized one-second bin data from each file
                        stimFileID = finalOptoStimFiles{jj,1};
                        stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
                        stimSpecField = 'cortical_LH';
                        for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
                            if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
                                stimCorticalS_Data_LH = AllSpecData.(stimSpecField).normS{kk,1};
                                F = AllSpecData.(stimSpecField).F{kk,1};
                                T = round(AllSpecData.(stimSpecField).T{kk,1},1);
                            end
                        end
                        stimStartTimeIndex = find(T == round(finalOptoStimStartTimes(jj,1),1));
                        stimStartTimeIndex = stimStartTimeIndex(1);
                        stimDurationIndex = find(T == round(finalOptoStimEndTimes(jj,1),1));
                        stimDurationIndex = stimDurationIndex(end);
                        stimCortS_Vals_LH = stimCorticalS_Data_LH(:,stimStartTimeIndex:stimDurationIndex);
                        stimCortZhold_LH = cat(3,stimCortZhold_LH,stimCortS_Vals_LH);
                    end
                    % cortical mean-subtract by first 3 seconds prior to stimulus
                    meanOptoStimCortS_LH = mean(stimCortZhold_LH,3);
                    baseOptoStimCortS_Vals_LH = mean(meanOptoStimCortS_LH(:,1:1.5*specSamplingRate),2);
                    baseMatrixOptoStimCortS_Vals_LH = baseOptoStimCortS_Vals_LH.*ones(size(meanOptoStimCortS_LH));
                    msOptoStimCortS_Vals_LH = (meanOptoStimCortS_LH - baseMatrixOptoStimCortS_Vals_LH);  
                    
                    % right hemispheres
                    stimCortZhold_RH = [];
                    for jj = 1:length(finalOptoStimFiles)
                        % load normalized one-second bin data from each file
                        stimFileID = finalOptoStimFiles{jj,1};
                        stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
                        stimSpecField = 'cortical_RH';
                        for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
                            if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
                                stimCorticalS_Data_RH = AllSpecData.(stimSpecField).normS{kk,1};
                                F = AllSpecData.(stimSpecField).F{kk,1};
                                T = round(AllSpecData.(stimSpecField).T{kk,1},1);
                            end
                        end
                        stimStartTimeIndex = find(T == round(finalOptoStimStartTimes(jj,1),1));
                        stimStartTimeIndex = stimStartTimeIndex(1);
                        stimDurationIndex = find(T == round(finalOptoStimEndTimes(jj,1),1));
                        stimDurationIndex = stimDurationIndex(end);
                        stimCortS_Vals_RH = stimCorticalS_Data_RH(:,stimStartTimeIndex:stimDurationIndex);
                        stimCortZhold_RH = cat(3,stimCortZhold_RH,stimCortS_Vals_RH);
                    end
                    % cortical mean-subtract by first 3 seconds prior to stimulus
                    meanOptoStimCortS_RH = mean(stimCortZhold_RH,3);
                    baseOptoStimCortS_Vals_RH = mean(meanOptoStimCortS_RH(:,1:1.5*specSamplingRate),2);
                    baseMatrixOptoStimCortS_Vals_RH = baseOptoStimCortS_Vals_RH.*ones(size(meanOptoStimCortS_RH));
                    msOptoStimCortS_Vals_RH = (meanOptoStimCortS_RH - baseMatrixOptoStimCortS_Vals_RH);  
        
                    %% save the data
                    T2 = -5:(1/specSamplingRate):15;
        
                    OptoStimData.P_NE.(solenoid).count = size(procOptoStimCortMUAData_LH,1);
            
                    OptoStimData.P_ACh.(solenoid).CBV.CBV = AChmeanOptoStimCBVData;
                    OptoStimData.P_ACh.(solenoid).CBV.CBVStD = AChstdOptoStimCBVData;
                    OptoStimData.P_NE.(solenoid).CBV.CBV = NEmeanOptoStimCBVData;
                    OptoStimData.P_NE.(solenoid).CBV.CBVStD = NEstdOptoStimCBVData;
                    OptoStimData.P_ACh.(solenoid).GFP.GFP= AChmeanOptoStimGFPData;
                    OptoStimData.P_ACh.(solenoid).GFP.GFPStD = AChstdOptoStimGFPData;
                    OptoStimData.P_NE.(solenoid).GFP.GFP= NEmeanOptoStimGFPData;
                    OptoStimData.P_NE.(solenoid).GFP.GFPStD = NEstdOptoStimGFPData;
        
                    OptoStimData.P_ACh.(solenoid).CBV.CBVRaw = AChprocOptoStimCBVData;
                    OptoStimData.P_NE.(solenoid).CBV.CBVRaw = NEprocOptoStimCBVData;
                    OptoStimData.P_ACh.(solenoid).GFP.GFPRaw = AChprocOptoStimGFPData;
                    OptoStimData.P_NE.(solenoid).GFP.GFPRaw = NEprocOptoStimGFPData;
            
                    OptoStimData.cortical_LH.(solenoid).MUA.corticalData = meanOptoStimCortMUAData_LH;
                    OptoStimData.cortical_LH.(solenoid).MUA.corticalStD = stdOptoStimCortMUAData_LH;
                    OptoStimData.cortical_LH.(solenoid).Gam.corticalData = meanOptoStimCortGamData_LH;
                    OptoStimData.cortical_LH.(solenoid).Gam.corticalStD = stdOptoStimCortGamData_LH;
        
                    OptoStimData.cortical_RH.(solenoid).MUA.corticalData = meanOptoStimCortMUAData_RH;
                    OptoStimData.cortical_RH.(solenoid).MUA.corticalStD = stdOptoStimCortMUAData_RH;
                    OptoStimData.cortical_RH.(solenoid).Gam.corticalData = meanOptoStimCortGamData_RH;
                    OptoStimData.cortical_RH.(solenoid).Gam.corticalStD = stdOptoStimCortGamData_RH;
        
                    OptoStimData.Pupil.(solenoid).Diameter.DiameterData = meanOptoStimPupilData;
                    OptoStimData.Pupil.(solenoid).Diameter.DiameterStD = stdOptoStimPupilData;
                    OptoStimData.Pupil.(solenoid).Diameter.DiameterRaw = procOptoStimPupilData;
        
                    OptoStimData.cortical_LH.(solenoid).timeVector = timeVector;
                    OptoStimData.cortical_LH.(solenoid).ECoG.corticalS_LH = msOptoStimCortS_Vals_LH;
                    OptoStimData.cortical_RH.(solenoid).ECoG.corticalS_RH = msOptoStimCortS_Vals_RH;
                    OptoStimData.cortical_LH.(solenoid).ECoG.T = T2;
                    OptoStimData.cortical_LH.(solenoid).ECoG.F = F;     
end
