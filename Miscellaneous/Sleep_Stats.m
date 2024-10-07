clear all
close all
clc


% find and load RestData.mat struct
    trainingDataFileStruct = dir('*_TrainingData.mat');
    trainingDataFile = {trainingDataFileStruct.name}';

for Tn = 1:1:size(trainingDataFile)
    trainDataFileID = char(trainingDataFile(Tn,1));
    load(trainDataFileID)
    ArousalState = trainingTable.behavState;
    Awake = zeros(size(ArousalState,1),1);
    NREM = zeros(size(ArousalState,1),1);
    REM = zeros(size(ArousalState,1),1);

        for SL = 1:1:size(ArousalState,1)
            if ArousalState(SL) == "Not Sleep"
                Awake(SL) = 1;
            elseif ArousalState(SL) == "NREM Sleep"
                NREM(SL) = 1;
            elseif ArousalState(SL) == "REM Sleep"
                REM(SL) = 1;
            end
        end
        if Tn == 1
            Awake_Total = Awake;
            NREM_Total = NREM;
            REM_Total = REM;
        else
            Awake_Total = cat(2,Awake_Total,Awake);
            NREM_Total = cat(2,NREM_Total,NREM);
            REM_Total = cat(2,REM_Total,REM);
        end
end

    Awake_Sum = sum(Awake_Total)/size(ArousalState,1);
    NREM_Sum = sum(NREM_Total)/size(ArousalState,1);
    REM_Sum = sum(REM_Total)/size(ArousalState,1);
    Varo = [Awake_Sum;NREM_Sum;REM_Sum];

    A_Mean = mean(Varo,2)';
    A_name = ["Awake","NREM","REM"];
    figure;
    p = piechart(A_Mean);
    p.Names = ["Awake","NREM","REM"];
  
        p.ColorOrder = [0 0 0;0 0 1;1 0 0];


Awake_std =  std(Awake_Sum);
NREM_std =  std(NREM_Sum);
REM_std =  std(REM_Sum);
%%

for Tn = 1:1:size(trainingDataFile)
    trainDataFileID = char(trainingDataFile(Tn,1));
    load(trainDataFileID)
    ArousalState = trainingTable.behavState;
    Awake = nan(size(ArousalState,1),1);
    NREM = nan(size(ArousalState,1),1);
    REM = nan(size(ArousalState,1),1);

        for SL = 1:1:size(ArousalState,1)
            if ArousalState(SL) == "Not Sleep"
                Awake(SL) = 1;
            elseif ArousalState(SL) == "NREM Sleep"
                NREM(SL) = 1;
            elseif ArousalState(SL) == "REM Sleep"
                REM(SL) = 1;
            end
        end
        if Tn == 1
            Awake_Total = Awake;
            NREM_Total = NREM;
            REM_Total = REM;
        else
            Awake_Total = cat(2,Awake_Total,Awake);
            NREM_Total = cat(2,NREM_Total,NREM);
            REM_Total = cat(2,REM_Total,REM);
        end
end    
AwakeYvals = 1 * Awake_Total; 
NREMYvals = 1 * NREM_Total; 
REMYvals = 1 * REM_Total; 

SleepDummy = 1:5:3120;


figure;
subplotn = 1;
for SJ = 1:1:14%size(AwakeYvals,2)
    subplot(14,1,subplotn);

    s1 = scatter(SleepDummy',AwakeYvals(:,SJ),'s','filled','MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on
    s2 = scatter(SleepDummy',NREMYvals(:,SJ),'s','MarkerFaceColor','b','MarkerEdgeColor','b');
    s3 = scatter(SleepDummy',REMYvals(:,SJ),'s','MarkerFaceColor','r','MarkerEdgeColor','r');
    xlim([1 3120]);
    ylim([0.99 1.01]);
    xticks([])
    set(gca,'TickLength',[0,0])
    yticks([])
    subplotn = subplotn+1;
end


figure;
subplotn = 1;
for SJ = 14+1:1:14+14%size(AwakeYvals,2)
    subplot(14,1,subplotn);

    s1 = scatter(SleepDummy',AwakeYvals(:,SJ),'s','filled','MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on
    s2 = scatter(SleepDummy',NREMYvals(:,SJ),'s','MarkerFaceColor','b','MarkerEdgeColor','b');
    s3 = scatter(SleepDummy',REMYvals(:,SJ),'s','MarkerFaceColor','r','MarkerEdgeColor','r');
    xlim([1 3120]);
    ylim([0.99 1.01]);
    xticks([])
    set(gca,'TickLength',[0,0])
    yticks([])
    subplotn = subplotn+1;
end

figure;
subplotn = 1;
for SJ = 28+1:1:28+14%size(AwakeYvals,2)
    subplot(14,1,subplotn);

    s1 = scatter(SleepDummy',AwakeYvals(:,SJ),'s','filled','MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on
    s2 = scatter(SleepDummy',NREMYvals(:,SJ),'s','MarkerFaceColor','b','MarkerEdgeColor','b');
    s3 = scatter(SleepDummy',REMYvals(:,SJ),'s','MarkerFaceColor','r','MarkerEdgeColor','r');
    xlim([1 3120]);
    ylim([0.99 1.01]);
    xticks([])
    set(gca,'TickLength',[0,0])
    yticks([])
    subplotn = subplotn+1;
end

figure;
subplotn = 1;
for SJ = 42+1:1:42+14%size(AwakeYvals,2)
    subplot(14,1,subplotn);

    s1 = scatter(SleepDummy',AwakeYvals(:,SJ),'s','filled','MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on
    s2 = scatter(SleepDummy',NREMYvals(:,SJ),'s','MarkerFaceColor','b','MarkerEdgeColor','b');
    s3 = scatter(SleepDummy',REMYvals(:,SJ),'s','MarkerFaceColor','r','MarkerEdgeColor','r');
    xlim([1 3120]);
    ylim([0.99 1.01]);
    xticks([])
    set(gca,'TickLength',[0,0])
    yticks([])
    subplotn = subplotn+1;
end