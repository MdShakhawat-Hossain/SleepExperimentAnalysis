% load('NEACh007_230719_14_52_21_ProcData.mat')
% function LocomotionTrigChange()
%% Run triggered averaging
RunDuration=[5 10 15 30 45];
RunLengths={'five_second_events','ten_second_events','fifteen_second_events','thirty_second_events','fortyfive_second_events'};
LeadTime=5;    
[~,~,new_T_run]=MotionTrigger_Fiber(ProcData.data.binForceSensor,ProcData.notes.dsFs,RunDuration(1),LeadTime);
EventLengths=(new_T_run(2,:)-new_T_run(1,:))/ProcData.notes.dsFs;
for durNum=1:length(RunDuration)
    LocomotionEvents=[];
    RunLabel=RunLengths{durNum};
    
    if durNum==length(RunDuration)
        EventFind=EventLengths>=RunDuration(durNum);
        EventInds=new_T_run(:,EventFind);
    else
        EventFind=EventLengths>=RunDuration(durNum) & EventLengths<RunDuration(durNum+1);
        EventInds=new_T_run(:,EventFind);
    end
    
    for k=1:size(EventInds,2)
        StartInd=EventInds(1,k);
        if durNum==length(RunDuration)
            EndInd=EventInds(1,k)+(ProcData.notes.dsFs*60);
        else
            EndInd=EventInds(1,k)+(ProcData.notes.dsFs*RunDuration(durNum+1));
        end
        % if StartInd>0
        %     if EndInd<length(ZscoredFiberData)+1
        %         Baseline=mean(ZscoredFiberData((StartInd:(StartInd+(2*ProcData.notes.dsFs))),:),1);
        %         BaselineRaw=mean(LowPassData((StartInd:(StartInd+(2*ProcData.notes.dsFs))),:),1);
        %         for j=1:size(Baseline,2)
        %             LocomotionEvents(:,j,k)=ZscoredFiberData((StartInd:EndInd),j)-Baseline(j);
        %             LocomotionEventsRaw(:,j,k)=LowPassData((StartInd:EndInd),j)-BaselineRaw(j);
        %         end
        %     end
        % end
    end
    
    % if ~isempty(LocomotionEvents)
    %     % AverageLocomotionResponse=mean(LocomotionEvents,3);
    %     % StanDevLocomotion=std(LocomotionEvents,0,3);
    %     % MedianLocomotionResponse=median(LocomotionEvents,3);
    %     % BaseResp=mean(AverageLocomotionResponse((1:(3*ProcData.notes.dsFs)),:),1);
    %     % PeakVal=max(AverageLocomotionResponse,[],1);
    %     % for chanNum=1:size(AverageLocomotionResponse,2)
    %     %     PeakTime(chanNum)=(find(AverageLocomotionResponse(:,chanNum)==PeakVal(chanNum))-ProcData.Params.StartPad)/ProcData.notes.dsFs;
    %     % end
    %     % PeakInds=((ProcData.Params.StartPad+RunDuration(durNum)*ProcData.notes.dsFs)-(3*ProcData.notes.dsFs)):(ProcData.Params.StartPad+RunDuration(durNum)*ProcData.notes.dsFs);
    %     % PeakResp=mean(AverageLocomotionResponse(PeakInds,:),1);
    % 
    %     % ProcData.WheelData.(RunLabel).RunInds=EventInds;
    %     % ProcData.WheelData.(RunLabel).RunLengths=(EventInds(2,:)-EventInds(1,:))/ProcData.notes.dsFs;
    %     % ProcData.AveragedData.(RunLabel).AvgResp=AverageLocomotionResponse;
    %     % ProcData.AveragedData.(RunLabel).StdResp=StanDevLocomotion;
    %     % ProcData.AveragedData.(RunLabel).MedResp=MedianLocomotionResponse;
    %     % ProcData.AveragedData.(RunLabel).BaseResp=BaseResp;
    %     % ProcData.AveragedData.(RunLabel).PeakResp=PeakResp;
    %     % ProcData.AveragedData.(RunLabel).PeakVal=PeakVal;
    %     % ProcData.AveragedData.(RunLabel).PeakTime=PeakTime;
    %     % ProcData.LocomotionEvokedData.(RunLabel).OpticalData=LocomotionEvents;
    %     % ProcData.LocomotionEvokedData.(RunLabel).OpticalDataRaw=LocomotionEventsRaw;
    % else
    %     % ProcData.WheelData.(RunLabel).RunInds=[];
    %     % ProcData.WheelData.(RunLabel).RunLengths=[];
    %     % ProcData.AveragedData.(RunLabel).AvgResp=[];
    %     % ProcData.AveragedData.(RunLabel).StdResp=[];
    %     % ProcData.AveragedData.(RunLabel).MedResp=[];
    %     % ProcData.AveragedData.(RunLabel).BaseResp=[];
    %     % ProcData.AveragedData.(RunLabel).PeakResp=[];
    %     % ProcData.AveragedData.(RunLabel).PeakVal=[];
    %     % ProcData.AveragedData.(RunLabel).PeakTime=[];
    %     % ProcData.LocomotionEvokedData.(RunLabel).OpticalData=[];
    % end
    
%     plotTime=((1:length(ProcData.AveragedData.(RunLabel).AvgResp))-ProcData.Params.StartPad)/ProcData.notes.dsFs;
%     figure;plot(plotTime,ProcData.AveragedData.(RunLabel).AvgResp(:,(2:3))); xlim([-5 RunDuration(durNum)]); legend({'nNosCre GCaMP6s','TRITC Blood Volume'});
%     title('Average locomotion evoked response'); xlabel('Time (sec)'); ylabel('Z score');
    clear LocomotionEvents LocomotionEventsRaw
end

% plotTime=((1:length(ProcData.AveragedData.MedResp))-ProcData.Params.StartPad)/ProcData.notes.dsFs;
% figure;plot(plotTime,ProcData.AveragedData.MedResp(:,(2:3))); xlim([-5 15]); legend({'nNosCre GCaMP6s','TRITC Blood Volume'});title('Median Locomotion evoked response');xlabel('Time (sec)');ylabel('Z score');

% figure;plot(ProcData.Xcorr.Lags,ProcData.Xcorr.CBV_GCaMP);xlim([-10 10]); title('Cross-correlation blood volume x GCaMP6s'); xlabel('Time (sec)'); ylabel('Normalized corr coeff');
% 
% runPlot(1:length(ProcData.RawFiberData.ZScoreFiberData(:,3)))=0;
% for runNum=1:length(ProcData.WheelData.five_second_events.RunInds)
%     runPlot(ProcData.WheelData.five_second_events.RunInds(1,runNum):ProcData.WheelData.five_second_events.RunInds(2,runNum))=1.1*max(ProcData.RawFiberData.ZScoreFiberData(:,3));
% end
% runPlot(runPlot==0)=NaN;
% plotTime=(1:length(ProcData.RawFiberData.ZScoreFiberData(:,3)))/ProcData.notes.dsFs;
% 
% figure;plot(plotTime,ProcData.RawFiberData.ZScoreFiberData(:,(2:3)));title('Signal intensity changes');xlabel('Time (sec)');ylabel('Z-score'); 
% xlim([0 plotTime(end)]);
% hold on; scatter(plotTime,runPlot,'k','filled'); legend({'Ca2+ GCaMP6s','TRITC Blood volume','Locomotion events'});