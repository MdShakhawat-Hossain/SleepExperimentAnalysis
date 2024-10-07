function [AnalysisResults] = AnalyzeTRITCSleepProbability(IOS_animalIDs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the probability of arousal-state classification based on hemodynamic [TRITC] changes (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
LH_allCatLabels = [];
RH_allCatLabels = [];
LH_TRITCallCatMeans = [];
RH_TRITCallCatMeans = [];
%% extract data from each animal's sleep scoring results
for aa = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,aa};
    dataLoc = [rootFolder '/' animalID '/CombinedImaging/'];
    cd(dataLoc)
    % add this animal's scoring labels with the other animal'ss
    scoringResults = 'Forest_ScoringResults.mat';
    load(scoringResults,'-mat')
    baselinesFile = [animalID '_RestingBaselines.mat'];
    load(baselinesFile)
    LH_allCatLabels = cat(1,LH_allCatLabels,ScoringResults.alllabels);
    RH_allCatLabels = cat(1,RH_allCatLabels,ScoringResults.alllabels);
    % take the mean of each 5 second bin
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    binSize = 5;   % seconds
    numBins = 180;   % 15 minutes with 5 sec bins
    samplingRate = 30;   % Hz
    samplesPerBin = binSize*samplingRate;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        load(procDataFileID,'-mat')
        [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(fileDate);
        for cc = 1:numBins
            if cc == 1
                LH_TRITCbinSamples = ProcData.data.TRITC.LH(1:samplesPerBin);
                RH_TRITCbinSamples = ProcData.data.TRITC.RH(1:samplesPerBin);
            else
                LH_TRITCbinSamples = ProcData.data.TRITC.LH((cc - 1)*samplesPerBin + 1:cc*samplesPerBin);
                RH_TRITCbinSamples = ProcData.data.TRITC.RH((cc - 1)*samplesPerBin + 1:cc*samplesPerBin);
            end
            LH_TRITCallCatMeans = cat(1,LH_TRITCallCatMeans,mean(LH_TRITCbinSamples));
            RH_TRITCallCatMeans = cat(1,RH_TRITCallCatMeans,mean(RH_TRITCbinSamples));
        end
    end
end
% concatenage LH/RH labels and mean values.
allCatLabels = cat(1,LH_allCatLabels,RH_allCatLabels);
TRITCallCatMeans = cat(1,LH_TRITCallCatMeans,RH_TRITCallCatMeans);
%% change arousal-state labels to numbers
zAwake = 1; zNrem = 1; zRem = 1;
for zz = 1:length(allCatLabels)
    if strcmp(allCatLabels{zz,1},'Not Sleep') == true
        awakeTRITC(zAwake,1) = TRITCallCatMeans(zz,1); %#ok<*NASGU>
        zAwake = zAwake + 1;
    elseif strcmp(allCatLabels{zz,1},'NREM Sleep') == true
        nremTRITC(zNrem,1) = TRITCallCatMeans(zz,1);
        zNrem = zNrem + 1;
    elseif strcmp(allCatLabels{zz,1},'REM Sleep') == true
        remTRITC(zRem,1) = TRITCallCatMeans(zz,1);
        zRem = zRem + 1;
    end
end
%% put each mean and scoring label into a cell
minTRITC = floor(min(TRITCallCatMeans));
maxTRITC = ceil(max(TRITCallCatMeans));
awakeBins = minTRITC:1:maxTRITC;
cutDown = abs(minTRITC - (-35));
cutUp = maxTRITC - 115;
probBinLabels = cell(length(minTRITC:1:maxTRITC),1);
probBinMeans = cell(length(minTRITC:1:maxTRITC),1);
discBins = discretize(TRITCallCatMeans,awakeBins);
for dd = 1:length(discBins)
    probBinLabels{discBins(dd),1} = cat(1,probBinLabels{discBins(dd),1},{allCatLabels(dd,1)});
    probBinMeans{discBins(dd),1} = cat(1,probBinMeans{discBins(dd),1},{TRITCallCatMeans(dd,1)});
end
% condense the left edges of the histogram bins to -35:1:120
cutDownLabels = [];
cutDownMeans = [];
for ee = 1:cutDown
    cutDownLabels = cat(1,cutDownLabels,probBinLabels{ee,1});
    cutDownMeans = cat(1,cutDownMeans,probBinMeans{ee,1});
end
% condense the right edges of the histogram to -35:1:120
cutUpLabels = [];
cutUpMeans = [];
for ff = 1:cutUp + 1
    cutUpLabels = cat(1,cutUpLabels,probBinLabels{end - ff,1});
    cutUpMeans = cat(1,cutUpMeans,probBinMeans{end - ff,1});
end
% reconstruct array of labels based on new edges
finCatLabels = cat(1,{cutDownLabels},probBinLabels(cutDown + 1:end - (cutUp + 2)),{cutUpLabels});
% strcmp the bins and if the bin is asleep (NREM/REM) set to 0, else set 1
for gg = 1:length(finCatLabels)
    for hh = 1:length(finCatLabels{gg,1})
        if strcmp(finCatLabels{gg,1}{hh,1}{1,1},'Not Sleep') == true
            awakeProbEvents{gg,1}(hh,1) = 1; %#ok<*AGROW>
        else
            awakeProbEvents{gg,1}(hh,1) = 0;
        end
    end
end
% strcmp the bins and if the bin is not in NREM (Awake/REM) set to 0, else set 1
for gg = 1:length(finCatLabels)
    for hh = 1:length(finCatLabels{gg,1})
        if strcmp(finCatLabels{gg,1}{hh,1}{1,1},'NREM Sleep') == true
            nremProbEvents{gg,1}(hh,1) = 1;
        else
            nremProbEvents{gg,1}(hh,1) = 0;
        end
    end
end
% strcmp the bins and if the bin is not in REM (Awake/NREM) set to 0, else set 1
for gg = 1:length(finCatLabels)
    for hh = 1:length(finCatLabels{gg,1})
        if strcmp(finCatLabels{gg,1}{hh,1}{1,1},'REM Sleep') == true
            remProbEvents{gg,1}(hh,1) = 1;
        else
            remProbEvents{gg,1}(hh,1) = 0;
        end
    end
end
% take probability of each bin
for ii = 1:length(awakeProbEvents)
    awakeProbPerc(ii,1) = sum(awakeProbEvents{ii,1})/length(awakeProbEvents{ii,1})*100;
    nremProbPerc(ii,1) = sum(nremProbEvents{ii,1})/length(nremProbEvents{ii,1})*100;
    remProbPerc(ii,1) = sum(remProbEvents{ii,1})/length(remProbEvents{ii,1})*100;
end
% save results
AnalysisResults.TRITCSleepProbability.TRITCCatMeans = TRITCallCatMeans;
AnalysisResults.TRITCSleepProbability.awakeProbPerc = awakeProbPerc;
AnalysisResults.TRITCSleepProbability.nremProbPerc = nremProbPerc;
AnalysisResults.TRITCSleepProbability.remProbPerc = remProbPerc;
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults')

end
