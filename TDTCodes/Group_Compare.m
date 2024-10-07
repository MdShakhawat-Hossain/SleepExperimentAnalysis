clear all; close all; clc

NE004 = load('C:\Users\vsn5034\Documents\Victoria\2022-06-30_GRAB-NE_Mouse_004_waterbath\FiberprocessedData.mat');
NE005 = load('C:\Users\vsn5034\Documents\Victoria\2022-06-30_GRAB-NE_Mouse_005_waterbath\FiberprocessedData.mat');

%%
   %% plot to compare the Gas Conditions during rest
    Air.restevent_465 = [NE004.Air.restevent_465,NE005.Air.restevent_465];
    O2.restevent_465 = [NE004.O2.restevent_465,NE005.O2.restevent_465];
   
   
    Air.restevent_560 = [NE004.Air.restevent_560,NE005.Air.restevent_560];
    O2.restevent_560 = [NE004.O2.restevent_560,NE005.O2.restevent_560];
   
   
    figure; 
    subplot(211)
    plot(ones(size(Air.restevent_465)),Air.restevent_465,'o');
    hold on
    plot(2*ones(size(O2.restevent_465)),O2.restevent_465,'o');
    boxplot([Air.restevent_465,O2.restevent_465]', [ones(size(Air.restevent_465)),2*ones(size(O2.restevent_465))])
    hold off
    ylabel('Average \DeltaF/F');
    xlim([0 3])
    xticks([1 2])
    xticklabels({'Air','O2'})
    title('Gas triggered change in GRAB-NE fluoresence during resting period');
        
    subplot(212)
    plot(ones(size(Air.restevent_560)),Air.restevent_560,'o');
    hold on
    plot(2*ones(size(O2.restevent_560)),O2.restevent_560,'o');
    boxplot([Air.restevent_560,O2.restevent_560]', [ones(size(Air.restevent_560)),2*ones(size(O2.restevent_560))])
    hold off
    xlabel('Gas');
    ylabel('Average \DeltaF/F');
    xlim([0 3])
    xticks([1 2])
    xticklabels({'Air','O2'})
    title('Event triggered change in CBV during resting period');  
    
   %% Anova
    Y_data_465 = [Air.restevent_465,O2.restevent_465];
    Grp_data_465 = [ones(size(Air.restevent_465)),2*ones(size(O2.restevent_465))];
    p_anova_465 = anova1(Y_data_465',Grp_data_465');
    
    Y_data_560 = [Air.restevent_560,O2.restevent_560];
    Grp_data_560 = [ones(size(Air.restevent_560)),2*ones(size(O2.restevent_560))];
    p_anova_560 = anova1(Y_data_560',Grp_data_560');
    
    
   %% plot to compare the Gas Conditions during rest in last 5 min
   Air.restevent_465_5min = [NE004.Air.restevent_465_5min,NE005.Air.restevent_465_5min];
   O2.restevent_465_5min = [NE004.O2.restevent_465_5min,NE005.O2.restevent_465_5min];
   
   
    Air.restevent_560_5min = [NE004.Air.restevent_560_5min,NE005.Air.restevent_560_5min];
   O2.restevent_560_5min = [NE004.O2.restevent_560_5min,NE005.O2.restevent_560_5min];
   
   
    figure; 
    subplot(211)
    plot(ones(size(Air.restevent_465_5min)),Air.restevent_465_5min,'o');
    hold on
    plot(2*ones(size(O2.restevent_465_5min)),O2.restevent_465_5min,'o');
    boxplot([Air.restevent_465_5min,O2.restevent_465_5min]', [ones(size(Air.restevent_465_5min)),2*ones(size(O2.restevent_465_5min))])
    hold off
    ylabel('Average \DeltaF/F');
    xlim([0 3])
    xticks([1 2])
    xticklabels({'Air','O2'})
    title('Gas triggered change in GRAB-NE fluoresence during resting period - 5 min');
        
    subplot(212)
    plot(ones(size(Air.restevent_560_5min)),Air.restevent_560_5min,'o');
    hold on
    plot(2*ones(size(O2.restevent_560_5min)),O2.restevent_560_5min,'o');
    boxplot([Air.restevent_560_5min,O2.restevent_560_5min]', [ones(size(Air.restevent_560_5min)),2*ones(size(O2.restevent_560_5min))])
    hold off
    xlabel('Gas');
    ylabel('Average \DeltaF/F');
    xlim([0 3])
    xticks([1 2])
    xticklabels({'Air','O2'})
    title('Event triggered change in CBV during resting period - 5 min');  
    
   %% Anova
    Y_data_465_5min = [Air.restevent_465_5min,O2.restevent_465_5min];
    Grp_data_465_5min = [ones(size(Air.restevent_465_5min)),2*ones(size(O2.restevent_465_5min))];
    p_anova_465_5min = anova1(Y_data_465_5min',Grp_data_465_5min');
    
    Y_data_560_5min = [Air.restevent_560_5min,O2.restevent_560_5min];
    Grp_data_560_5min = [ones(size(Air.restevent_560_5min)),2*ones(size(O2.restevent_560_5min))];
    p_anova_560_5min = anova1(Y_data_560_5min',Grp_data_560_5min');
    
    
  %% Find movement and separate those from run
  
    speedaccBinary = speedProcess_Doric(AnalogWheel,Params.DataFs);
    
    Speed_Out = findSpeedSeg(speedaccBinary, Params.DataFs,30, 0);
    
    % mean value for all resting period
    for restidx = 1:1:length(Speed_Out.Rest(:,1))
        restevent_465(restidx) = mean(FP465.dF(Speed_Out.Rest(restidx,1):Speed_Out.Rest(restidx,2))); 
    end
    average_465.rest = mean(restevent_465);

    for runidx = 1:1:length(Speed_Out.Run(:,1))
        runevent_465(runidx) = mean(FP465.dF(Speed_Out.Run(runidx,1):Speed_Out.Run(runidx,2))); 
    end
    average_465.run = mean(runevent_465);
    %      
    for restidx = 1:1:length(Speed_Out.Rest(:,1))
        restevent_560(restidx) = mean(FP560.dF(Speed_Out.Rest(restidx,1):Speed_Out.Rest(restidx,2))); 
    end
    average_560.rest = mean(restevent_560);
            
    for runidx = 1:1:length(Speed_Out.Run(:,1))
        runevent_560(runidx) = mean(FP560.dF(Speed_Out.Run(runidx,1):Speed_Out.Run(runidx,2))); 
    end
    average_560.run = mean(runevent_560);
    % plot the event triggered changes in fluoresence
    
    figure; 
    subplot(211)
    plot(ones(size(restevent_465)),restevent_465,'o');
    hold on
    plot(2*ones(size(runevent_465)),runevent_465,'o');
    boxplot([restevent_465;runevent_465]', {'Rest','Run'})
    hold off 
    ylabel('Average \DeltaF/F');
    title('Event triggered change in GRAB-NE fluoresence');
        
    subplot(212)
    plot(ones(size(restevent_560)),restevent_560,'o');
    hold on
    plot(2*ones(size(runevent_560)),runevent_560,'o');
    boxplot([restevent_560;runevent_560]', {'Rest','Run'})
    hold off
    xlabel('Event');
    ylabel('Average \DeltaF/F');
    title('Event triggered change in CBV');
  
    %% only the events that happened during gas conditions
    
    Xvalues = Xvalues.*Params.DataFs*60;
    Vals_Air = [1,2,5,6,9,10];
    Vals_O2 = [3,4,7,8,11,12];
    arsIdx=1;
    arnIdx=1;
    orsIdx=1;
    ornIdx=1;
    for kIdx = 1:1:length(Speed_Out.Rest(:,1))
        for gasevent = 1:1:length(Vals_Air)/2
            if Speed_Out.Rest(kIdx,1) >= Xvalues(Vals_Air((gasevent*2)-1)) && Speed_Out.Rest(kIdx,2) <= Xvalues(Vals_Air(gasevent*2))
                RestIdx.Air(arsIdx,:) = Speed_Out.Rest(kIdx,:);
                arsIdx = arsIdx+1;
            end
            if Speed_Out.Run(kIdx,1) >= Xvalues(Vals_Air((gasevent*2)-1)) && Speed_Out.Rest(kIdx,2) <= Xvalues(Vals_Air(gasevent*2))
                RunIdx.Air(arnIdx,:) = Speed_Out.Run(kIdx,:);
                arnIdx = arnIdx+1;
            end
        
            if Speed_Out.Rest(kIdx,1) >= Xvalues(Vals_O2((gasevent*2)-1)) && Speed_Out.Rest(kIdx,2) <= Xvalues(Vals_O2(gasevent*2))
                RestIdx.O2(orsIdx,:) = Speed_Out.Rest(kIdx,:);
                orsIdx = orsIdx+1;
            end
            if Speed_Out.Run(kIdx,1) >= Xvalues(Vals_O2((gasevent*2)-1)) && Speed_Out.Rest(kIdx,2) <= Xvalues(Vals_O2(gasevent*2))
                RunIdx.O2(ornIdx,:) = Speed_Out.Run(kIdx,:);
                ornIdx = ornIdx+1;
            end
        end
    end
     %% gas event triggered average
     %%%  Air ----------
     % mean value for resting period during Air
     
    Air.restevent_465 = [NE004.Air.restevent_465,NE005.Air.restevent_465];
    Air.restevent_560 = [NE004.Air.restevent_560,NE005.Air.restevent_560];
    
    for restidx = 1:1:length(RestIdx.Air(:,1))
        Air.restevent_465(restidx) = mean(FP465.dF(RestIdx.Air(restidx,1):RestIdx.Air(restidx,2))); 
    end
    Air.average_465.rest = mean(Air.restevent_465);
    
    for restidx = 1:1:length(RestIdx.Air(:,1))
        Air.restevent_560(restidx) = mean(FP560.dF(RestIdx.Air(restidx,1):RestIdx.Air(restidx,2))); 
    end
    Air.average_560.rest = mean(Air.restevent_560);
    % mean value for running period during Air
    for runidx = 1:1:length(RunIdx.Air(:,1))
        Air.runevent_465(runidx) = mean(FP465.dF(RunIdx.Air(runidx,1):RunIdx.Air(runidx,2))); 
    end
    Air.average_465.run = mean(Air.runevent_465);
    
    for runidx = 1:1:length(RunIdx.Air(:,1))
        Air.runevent_560(runidx) = mean(FP560.dF(RunIdx.Air(runidx,1):RunIdx.Air(runidx,2))); 
    end
    Air.average_560.run = mean(Air.runevent_560);
    % plot the event triggered changes in fluoresence
    
    figure; 
    subplot(211)
    plot(ones(size(Air.restevent_465)),Air.restevent_465,'o');
    hold on
    plot(2*ones(size(Air.runevent_465)),Air.runevent_465,'o');
    boxplot([Air.restevent_465,Air.runevent_465]',[ones(size(Air.restevent_465)),2*ones(size(Air.runevent_465))]) 
    hold off 
    ylabel('Average \DeltaF/F');
    xlim([0 3])
    xticks([1 2])
    xticklabels({'Rest','Run'})
    title('Event triggered change in GRAB-NE fluoresence for Air');
        
    subplot(212)
    plot(ones(size(Air.restevent_560)),Air.restevent_560,'o');
    hold on
    plot(2*ones(size(Air.runevent_560)),Air.runevent_560,'o');
    boxplot([Air.restevent_560,Air.runevent_560]', [ones(size(Air.restevent_560)),2*ones(size(Air.runevent_560))]) 
    hold off
    xlabel('Event');
    ylabel('Average \DeltaF/F');
    xlim([0 3])
    xticks([1 2])
    xticklabels({'Rest','Run'})
    title('Event triggered change in CBV for Air');
    
    %%%  O2 ----------
     % mean value for resting period during O2
    
    O2.restevent_465 = [NE004.O2.restevent_465,NE005.O2.restevent_465];
    O2.restevent_560 = [NE004.O2.restevent_560,NE005.O2.restevent_560];
     
    for restidx = 1:1:length(RestIdx.O2(:,1))
        O2.restevent_465(restidx) = mean(FP465.dF(RestIdx.O2(restidx,1):RestIdx.O2(restidx,2))); 
    end
    O2.average_465.rest = mean(O2.restevent_465);
    
    for restidx = 1:1:length(RestIdx.O2(:,1))
        O2.restevent_560(restidx) = mean(FP560.dF(RestIdx.O2(restidx,1):RestIdx.O2(restidx,2))); 
    end
    O2.average_560.rest = mean(O2.restevent_560);
    % mean value for running period during O2
    for runidx = 1:1:length(RunIdx.O2(:,1))
        O2.runevent_465(runidx) = mean(FP465.dF(RunIdx.O2(runidx,1):RunIdx.O2(runidx,2))); 
    end
    O2.average_465.run = mean(O2.runevent_465);
    
    for runidx = 1:1:length(RunIdx.O2(:,1))
        O2.runevent_560(runidx) = mean(FP560.dF(RunIdx.O2(runidx,1):RunIdx.O2(runidx,2))); 
    end
    O2.average_560.run = mean(O2.runevent_560);
    % plot the event triggered changes in fluoresence
    
    figure; 
    subplot(211)
    plot(ones(size(O2.restevent_465)),O2.restevent_465,'o');
    hold on
    plot(2*ones(size(O2.runevent_465)),O2.runevent_465,'o');
    boxplot([O2.restevent_465,O2.runevent_465]', [ones(size(O2.restevent_465)),2*ones(size(O2.runevent_465))]) 
    hold off 
    ylabel('Average \DeltaF/F');
    xlim([0 3])
    xticks([1 2])
    xticklabels({'Rest','Run'})
    title('Event triggered change in GRAB-NE fluoresence for O2');
        
    subplot(212)
    plot(ones(size(O2.restevent_560)),O2.restevent_560,'o');
    hold on
    plot(2*ones(size(O2.runevent_560)),O2.runevent_560,'o');
    boxplot([O2.restevent_560,O2.runevent_560]', [ones(size(O2.restevent_560)),2*ones(size(O2.runevent_560))]) 
    hold off 
    ylabel('Average \DeltaF/F');
    xlabel('Event')
    xlim([0 3])
    xticks([1 2])
    xticklabels({'Rest','Run'})
    title('Event triggered change in CBV for O2');
    
    
    %% Anova
    Y_data_560 = [Air.restevent_560,O2.restevent_560];
    Grp_data_560 = [ones(size(Air.restevent_560)),2*ones(size(O2.restevent_560))];
    p_anova_560 = anova1(Y_data_560',Grp_data_560');
   
    Y_data_465 = [Air.restevent_465,O2.restevent_465];
    Grp_data_465 = [ones(size(Air.restevent_465)),2*ones(size(O2.restevent_465))];
    p_anova_465 = anova1(Y_data_465',Grp_data_465');  
    
    
    
    