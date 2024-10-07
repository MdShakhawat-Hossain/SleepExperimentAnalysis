function [AnalysisResults] = Fig7(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 7 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
% colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
colorRest = [(0/256),(166/256),(81/256)];
colorWhisk = [(31/256),(120/256),(179/256)];
% colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
colorAlert = [(255/256),(191/256),(0/256)];
colorAsleep = [(0/256),(128/256),(255/256)];
colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
IOS_animalIDs = {'T281','T282','T284','T285'};
behavFields = {'Rest','NREM','REM','Awake','All'};
% behavFields2 = {'Rest','NREM','REM','Awake','All'};
behavFields3 = {'Rest','Whisk','NREM','REM','Awake','All'};
dataTypes = {'HbT','gammaBandPower'};
%% average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.Coherr = [];
for a = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        % create the behavior folder for the first iteration of the loop
        if isfield(data.Coherr,behavField) == false
            data.Coherr.(behavField) = [];
        end
        for c = 1:length(dataTypes)
            dataType = dataTypes{1,c};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).Coherence.(behavField).(dataType).C) == false
                % create the data type folder for the first iteration of the loop
                if isfield(data.Coherr.(behavField),dataType) == false
                    data.Coherr.(behavField).(dataType).C = [];
                    data.Coherr.(behavField).(dataType).f = [];
                    data.Coherr.(behavField).(dataType).confC = [];
                end
                % concatenate C/f for existing data - exclude any empty sets
                data.Coherr.(behavField).(dataType).C = cat(2,data.Coherr.(behavField).(dataType).C,AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                data.Coherr.(behavField).(dataType).f = cat(1,data.Coherr.(behavField).(dataType).f,AnalysisResults.(animalID).Coherence.(behavField).(dataType).f);
                data.Coherr.(behavField).(dataType).confC = cat(1,data.Coherr.(behavField).(dataType).confC,AnalysisResults.(animalID).Coherence.(behavField).(dataType).confC);
            end
        end
    end
end
% take mean/StD of C/f and determine confC line
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
        data.Coherr.(behavField).(dataType).meanC = mean(data.Coherr.(behavField).(dataType).C,2);
        data.Coherr.(behavField).(dataType).stdC = std(data.Coherr.(behavField).(dataType).C,0,2);
        data.Coherr.(behavField).(dataType).meanf = mean(data.Coherr.(behavField).(dataType).f,1);
        data.Coherr.(behavField).(dataType).maxConfC = geomean(data.Coherr.(behavField).(dataType).confC);
        data.Coherr.(behavField).(dataType).maxConfC_Y = ones(length(data.Coherr.(behavField).(dataType).meanf),1)*data.Coherr.(behavField).(dataType).maxConfC;
    end
end
%% power spectra during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
% for a = 1:length(IOS_animalIDs)
%     animalID = IOS_animalIDs{1,a};
%     for b = 1:length(behavFields)
%         behavField = behavFields{1,b};
%         for c = 1:length(dataTypes)
%             dataType = dataTypes{1,c};
%             data.PowerSpec.(behavField).(dataType).LH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).LH.S;
%             data.PowerSpec.(behavField).(dataType).LH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).LH.f;
%             data.PowerSpec.(behavField).(dataType).RH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).RH.S;
%             data.PowerSpec.(behavField).(dataType).RH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).RH.f;
%         end
%     end
% end
% find the peak of the resting PSD for each animal/hemisphere
% for a = 1:length(IOS_animalIDs)
%     for c = 1:length(dataTypes)
%         dataType = dataTypes{1,c};
%         data.PowerSpec.baseline.(dataType).LH{a,1} = max(data.PowerSpec.Rest.(dataType).LH.S{a,1});
%         data.PowerSpec.baseline.(dataType).RH{a,1} = max(data.PowerSpec.Rest.(dataType).RH.S{a,1});
%     end
% end
% DC-shift each animal/hemisphere/behavior PSD with respect to the resting peak
% for a = 1:length(IOS_animalIDs)
%     for dd = 1:length(behavFields)
%         behavField = behavFields{1,dd};
%         for j = 1:length(dataTypes)
%             dataType = dataTypes{1,j};
%             for ee = 1:size(data.PowerSpec.(behavField).(dataType).LH.S,2)
%                 data.PowerSpec.(behavField).(dataType).normLH{a,1} = (data.PowerSpec.(behavField).(dataType).LH.S{a,1})*(1/(data.PowerSpec.baseline.(dataType).LH{a,1}));
%                 data.PowerSpec.(behavField).(dataType).normRH{a,1} = (data.PowerSpec.(behavField).(dataType).RH.S{a,1})*(1/(data.PowerSpec.baseline.(dataType).RH{a,1}));
%             end
%         end
%     end
% end
% % concatenate the data from the left and right hemispheres - removes any empty data
% for ee = 1:length(behavFields)
%     behavField = behavFields{1,ee};
%     for ff = 1:length(dataTypes)
%         dataType = dataTypes{1,ff};
%         data.PowerSpec.(behavField).(dataType).cat_S = [];
%         data.PowerSpec.(behavField).(dataType).cat_f = [];
%         for z = 1:length(data.PowerSpec.(behavField).(dataType).normLH)
%             data.PowerSpec.(behavField).(dataType).cat_S = cat(2,data.PowerSpec.(behavField).(dataType).cat_S,data.PowerSpec.(behavField).(dataType).normLH{z,1},data.PowerSpec.(behavField).(dataType).normRH{z,1});
%             data.PowerSpec.(behavField).(dataType).cat_f = cat(1,data.PowerSpec.(behavField).(dataType).cat_f,data.PowerSpec.(behavField).(dataType).LH.f{z,1},data.PowerSpec.(behavField).(dataType).RH.f{z,1});
%         end
%     end
% end
% take mean/StD of S/f
% for h = 1:length(behavFields)
%     behavField = behavFields{1,h};
%     for j = 1:length(dataTypes)
%         dataType = dataTypes{1,j};
%         data.PowerSpec.(behavField).(dataType).meanCortS = mean(data.PowerSpec.(behavField).(dataType).cat_S,2);
%         data.PowerSpec.(behavField).(dataType).stdCortS = std(data.PowerSpec.(behavField).(dataType).cat_S,0,2);
%         data.PowerSpec.(behavField).(dataType).meanCortf = mean(data.PowerSpec.(behavField).(dataType).cat_f,1);
%     end
% end
%% Pearson's correlations during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.CorrCoef = [];
for a = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,a};
    for b = 1:length(behavFields3)
        behavField = behavFields3{1,b};
        % create the behavior folder for the first iteration of the loop
        if isfield(data.CorrCoef,behavField) == false
            data.CorrCoef.(behavField) = [];
        end
        for c = 1:length(dataTypes)
            dataType = dataTypes{1,c};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).meanR) == false
                % create the data type folder for the first iteration of the loop
                if isfield(data.CorrCoef.(behavField),dataType) == false
                    data.CorrCoef.(behavField).(dataType).meanRs = [];
                    data.CorrCoef.(behavField).(dataType).animalID = {};
                    data.CorrCoef.(behavField).(dataType).behavior = {};
                end
                % concatenate mean R and the animalID/behavior for statistics table
                data.CorrCoef.(behavField).(dataType).meanRs = cat(1,data.CorrCoef.(behavField).(dataType).meanRs,AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).meanR);
                data.CorrCoef.(behavField).(dataType).animalID = cat(1,data.CorrCoef.(behavField).(dataType).animalID,animalID);
                data.CorrCoef.(behavField).(dataType).behavior = cat(1,data.CorrCoef.(behavField).(dataType).behavior,behavField);
            end
        end
    end
end
% take mean/STD of R
for ee = 1:length(behavFields3)
    behavField = behavFields3{1,ee};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
        data.CorrCoef.(behavField).(dataType).meanR = mean(data.CorrCoef.(behavField).(dataType).meanRs,1);
        data.CorrCoef.(behavField).(dataType).stdR = std(data.CorrCoef.(behavField).(dataType).meanRs,0,1);
    end
end
%% statistics - generalized linear mixed-effects model
% HbT
% HbTtableSize = cat(1,data.CorrCoef.Rest.HbT.meanRs,data.CorrCoef.Whisk.HbT.meanRs,data.CorrCoef.NREM.HbT.meanRs,data.CorrCoef.REM.HbT.meanRs,...
%     data.CorrCoef.Awake.HbT.meanRs,data.CorrCoef.Sleep.HbT.meanRs,data.CorrCoef.All.HbT.meanRs);
% HbTTable = table('Size',[size(HbTtableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
% HbTTable.Mouse = cat(1,data.CorrCoef.Rest.HbT.animalID,data.CorrCoef.Whisk.HbT.animalID,data.CorrCoef.NREM.HbT.animalID,data.CorrCoef.REM.HbT.animalID,...
%     data.CorrCoef.Awake.HbT.animalID,data.CorrCoef.Sleep.HbT.animalID,data.CorrCoef.All.HbT.animalID);
% HbTTable.CorrCoef = cat(1,data.CorrCoef.Rest.HbT.meanRs,data.CorrCoef.Whisk.HbT.meanRs,data.CorrCoef.NREM.HbT.meanRs,data.CorrCoef.REM.HbT.meanRs,...
%     data.CorrCoef.Awake.HbT.meanRs,data.CorrCoef.Sleep.HbT.meanRs,data.CorrCoef.All.HbT.meanRs);
% HbTTable.Behavior = cat(1,data.CorrCoef.Rest.HbT.behavior,data.CorrCoef.Whisk.HbT.behavior,data.CorrCoef.NREM.HbT.behavior,data.CorrCoef.REM.HbT.behavior,...
%     data.CorrCoef.Awake.HbT.behavior,data.CorrCoef.Sleep.HbT.behavior,data.CorrCoef.All.HbT.behavior);
% HbTFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
% HbTStats = fitglme(HbTTable,HbTFitFormula);
% % gamma-band power
% gammatableSize = cat(1,data.CorrCoef.Rest.gammaBandPower.meanRs,data.CorrCoef.Whisk.gammaBandPower.meanRs,data.CorrCoef.NREM.gammaBandPower.meanRs,data.CorrCoef.REM.gammaBandPower.meanRs,...
%     data.CorrCoef.Awake.gammaBandPower.meanRs,data.CorrCoef.Sleep.gammaBandPower.meanRs,data.CorrCoef.All.gammaBandPower.meanRs);
% gammaTable = table('Size',[size(gammatableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
% gammaTable.Mouse = cat(1,data.CorrCoef.Rest.gammaBandPower.animalID,data.CorrCoef.Whisk.gammaBandPower.animalID,data.CorrCoef.NREM.gammaBandPower.animalID,data.CorrCoef.REM.gammaBandPower.animalID,...
%     data.CorrCoef.Awake.gammaBandPower.animalID,data.CorrCoef.Sleep.gammaBandPower.animalID,data.CorrCoef.All.gammaBandPower.animalID);
% gammaTable.CorrCoef = cat(1,data.CorrCoef.Rest.gammaBandPower.meanRs,data.CorrCoef.Whisk.gammaBandPower.meanRs,data.CorrCoef.NREM.gammaBandPower.meanRs,data.CorrCoef.REM.gammaBandPower.meanRs,...
%     data.CorrCoef.Awake.gammaBandPower.meanRs,data.CorrCoef.Sleep.gammaBandPower.meanRs,data.CorrCoef.All.gammaBandPower.meanRs);
% gammaTable.Behavior = cat(1,data.CorrCoef.Rest.gammaBandPower.behavior,data.CorrCoef.Whisk.gammaBandPower.behavior,data.CorrCoef.NREM.gammaBandPower.behavior,data.CorrCoef.REM.gammaBandPower.behavior,...
%     data.CorrCoef.Awake.gammaBandPower.behavior,data.CorrCoef.Sleep.gammaBandPower.behavior,data.CorrCoef.All.gammaBandPower.behavior);
% gammaFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
% gammaStats = fitglme(gammaTable,gammaFitFormula);
%% Fig. 7
summaryFigure = figure('Name','Fig7 (a-g)');
sgtitle('Figure 7 - Turner et al. 2020')
CC_xInds = ones(1,length(IOS_animalIDs));
CC_xInds2 = ones(1,length(data.CorrCoef.Awake.HbT.animalID));
% CC_xInds3 = ones(1,length(data.CorrCoef.Sleep.HbT.animalID));
%% average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.Coherr = [];
for a = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        % create the behavior folder for the first iteration of the loop
        if isfield(data.Coherr,behavField) == false
            data.Coherr.(behavField) = [];
        end
        for c = 1:length(dataTypes)
            dataType = dataTypes{1,c};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).Coherence.(behavField).(dataType).C) == false
                % create the data type folder for the first iteration of the loop
                if isfield(data.Coherr.(behavField),dataType) == false
                    data.Coherr.(behavField).(dataType).C = [];
                    data.Coherr.(behavField).(dataType).f = [];
                    data.Coherr.(behavField).(dataType).confC = [];
                end
                % concatenate C/f for existing data - exclude any empty sets
                data.Coherr.(behavField).(dataType).C = cat(2,data.Coherr.(behavField).(dataType).C,AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                data.Coherr.(behavField).(dataType).f = cat(1,data.Coherr.(behavField).(dataType).f,AnalysisResults.(animalID).Coherence.(behavField).(dataType).f);
                data.Coherr.(behavField).(dataType).confC = cat(1,data.Coherr.(behavField).(dataType).confC,AnalysisResults.(animalID).Coherence.(behavField).(dataType).confC);
            end
        end
    end
end
% take mean/StD of C/f and determine confC line
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
        data.Coherr.(behavField).(dataType).meanC = mean(data.Coherr.(behavField).(dataType).C,2);
        data.Coherr.(behavField).(dataType).stdC = std(data.Coherr.(behavField).(dataType).C,0,2);
        data.Coherr.(behavField).(dataType).meanf = mean(data.Coherr.(behavField).(dataType).f,1);
        data.Coherr.(behavField).(dataType).maxConfC = geomean(data.Coherr.(behavField).(dataType).confC);
        data.Coherr.(behavField).(dataType).maxConfC_Y = ones(length(data.Coherr.(behavField).(dataType).meanf),1)*data.Coherr.(behavField).(dataType).maxConfC;
    end
end
%% power spectra during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        for c = 1:length(dataTypes)
            dataType = dataTypes{1,c};
            data.PowerSpec.(behavField).(dataType).LH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).LH.S;
            data.PowerSpec.(behavField).(dataType).LH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).LH.f;
            data.PowerSpec.(behavField).(dataType).RH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).RH.S;
            data.PowerSpec.(behavField).(dataType).RH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).RH.f;
        end
    end
end
% find the peak of the resting PSD for each animal/hemisphere
for a = 1:length(IOS_animalIDs)
    for c = 1:length(dataTypes)
        dataType = dataTypes{1,c};
        data.PowerSpec.baseline.(dataType).LH{a,1} = max(data.PowerSpec.Rest.(dataType).LH.S{a,1});
        data.PowerSpec.baseline.(dataType).RH{a,1} = max(data.PowerSpec.Rest.(dataType).RH.S{a,1});
    end
end
% DC-shift each animal/hemisphere/behavior PSD with respect to the resting peak
for a = 1:length(IOS_animalIDs)
    for dd = 1:length(behavFields)
        behavField = behavFields{1,dd};
        for j = 1:length(dataTypes)
            dataType = dataTypes{1,j};
            for ee = 1:size(data.PowerSpec.(behavField).(dataType).LH.S,2)
                data.PowerSpec.(behavField).(dataType).normLH{a,1} = (data.PowerSpec.(behavField).(dataType).LH.S{a,1})*(1/(data.PowerSpec.baseline.(dataType).LH{a,1}));
                data.PowerSpec.(behavField).(dataType).normRH{a,1} = (data.PowerSpec.(behavField).(dataType).RH.S{a,1})*(1/(data.PowerSpec.baseline.(dataType).RH{a,1}));
            end
        end
    end
end
% % concatenate the data from the left and right hemispheres - removes any empty data
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
        data.PowerSpec.(behavField).(dataType).cat_S = [];
        data.PowerSpec.(behavField).(dataType).cat_f = [];
        for z = 1:length(data.PowerSpec.(behavField).(dataType).normLH)
            data.PowerSpec.(behavField).(dataType).cat_S = cat(2,data.PowerSpec.(behavField).(dataType).cat_S,data.PowerSpec.(behavField).(dataType).normLH{z,1},data.PowerSpec.(behavField).(dataType).normRH{z,1});
            data.PowerSpec.(behavField).(dataType).cat_f = cat(1,data.PowerSpec.(behavField).(dataType).cat_f,data.PowerSpec.(behavField).(dataType).LH.f{z,1},data.PowerSpec.(behavField).(dataType).RH.f{z,1});
        end
    end
end
% take mean/StD of S/f
for h = 1:length(behavFields)
    behavField = behavFields{1,h};
    for j = 1:length(dataTypes)
        dataType = dataTypes{1,j};
        data.PowerSpec.(behavField).(dataType).meanCortS = mean(data.PowerSpec.(behavField).(dataType).cat_S,2);
        data.PowerSpec.(behavField).(dataType).stdCortS = std(data.PowerSpec.(behavField).(dataType).cat_S,0,2);
        data.PowerSpec.(behavField).(dataType).meanCortf = mean(data.PowerSpec.(behavField).(dataType).cat_f,1);
    end
end
%% Pearson's correlations during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.CorrCoef = [];
for a = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,a};
    for b = 1:length(behavFields3)
        behavField = behavFields3{1,b};
        % create the behavior folder for the first iteration of the loop
        if isfield(data.CorrCoef,behavField) == false
            data.CorrCoef.(behavField) = [];
        end
        for c = 1:length(dataTypes)
            dataType = dataTypes{1,c};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).meanR) == false
                % create the data type folder for the first iteration of the loop
                if isfield(data.CorrCoef.(behavField),dataType) == false
                    data.CorrCoef.(behavField).(dataType).meanRs = [];
                    data.CorrCoef.(behavField).(dataType).animalID = {};
                    data.CorrCoef.(behavField).(dataType).behavior = {};
                end
                % concatenate mean R and the animalID/behavior for statistics table
                data.CorrCoef.(behavField).(dataType).meanRs = cat(1,data.CorrCoef.(behavField).(dataType).meanRs,AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).meanR);
                data.CorrCoef.(behavField).(dataType).animalID = cat(1,data.CorrCoef.(behavField).(dataType).animalID,animalID);
                data.CorrCoef.(behavField).(dataType).behavior = cat(1,data.CorrCoef.(behavField).(dataType).behavior,behavField);
            end
        end
    end
end
% take mean/STD of R
for ee = 1:length(behavFields3)
    behavField = behavFields3{1,ee};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
        data.CorrCoef.(behavField).(dataType).meanR = mean(data.CorrCoef.(behavField).(dataType).meanRs,1);
        data.CorrCoef.(behavField).(dataType).stdR = std(data.CorrCoef.(behavField).(dataType).meanRs,0,1);
    end
end
%% statistics - generalized linear mixed-effects model
% GCaMP7s
% GCaMP7stableSize = cat(1,data.CorrCoef.Rest.GCaMP7s.meanRs,data.CorrCoef.Whisk.GCaMP7s.meanRs,data.CorrCoef.NREM.GCaMP7s.meanRs,data.CorrCoef.REM.GCaMP7s.meanRs,...
%     data.CorrCoef.Awake.GCaMP7s.meanRs,data.CorrCoef.Sleep.GCaMP7s.meanRs,data.CorrCoef.All.GCaMP7s.meanRs);
% GCaMP7sTable = table('Size',[size(GCaMP7stableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
% GCaMP7sTable.Mouse = cat(1,data.CorrCoef.Rest.GCaMP7s.animalID,data.CorrCoef.Whisk.GCaMP7s.animalID,data.CorrCoef.NREM.GCaMP7s.animalID,data.CorrCoef.REM.GCaMP7s.animalID,...
%     data.CorrCoef.Awake.GCaMP7s.animalID,data.CorrCoef.Sleep.GCaMP7s.animalID,data.CorrCoef.All.GCaMP7s.animalID);
% GCaMP7sTable.CorrCoef = cat(1,data.CorrCoef.Rest.GCaMP7s.meanRs,data.CorrCoef.Whisk.GCaMP7s.meanRs,data.CorrCoef.NREM.GCaMP7s.meanRs,data.CorrCoef.REM.GCaMP7s.meanRs,...
%     data.CorrCoef.Awake.GCaMP7s.meanRs,data.CorrCoef.Sleep.GCaMP7s.meanRs,data.CorrCoef.All.GCaMP7s.meanRs);
% GCaMP7sTable.Behavior = cat(1,data.CorrCoef.Rest.GCaMP7s.behavior,data.CorrCoef.Whisk.GCaMP7s.behavior,data.CorrCoef.NREM.GCaMP7s.behavior,data.CorrCoef.REM.GCaMP7s.behavior,...
%     data.CorrCoef.Awake.GCaMP7s.behavior,data.CorrCoef.Sleep.GCaMP7s.behavior,data.CorrCoef.All.GCaMP7s.behavior);
% GCaMP7sFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
% GCaMP7sStats = fitglme(GCaMP7sTable,GCaMP7sFitFormula);
% % gamma-band power
% gammatableSize = cat(1,data.CorrCoef.Rest.GCaMP7s.meanRs,data.CorrCoef.Whisk.GCaMP7s.meanRs,data.CorrCoef.NREM.GCaMP7s.meanRs,data.CorrCoef.REM.GCaMP7s.meanRs,...
%     data.CorrCoef.Awake.GCaMP7s.meanRs,data.CorrCoef.Sleep.GCaMP7s.meanRs,data.CorrCoef.All.GCaMP7s.meanRs);
% gammaTable = table('Size',[size(gammatableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
% gammaTable.Mouse = cat(1,data.CorrCoef.Rest.GCaMP7s.animalID,data.CorrCoef.Whisk.GCaMP7s.animalID,data.CorrCoef.NREM.GCaMP7s.animalID,data.CorrCoef.REM.GCaMP7s.animalID,...
%     data.CorrCoef.Awake.GCaMP7s.animalID,data.CorrCoef.Sleep.GCaMP7s.animalID,data.CorrCoef.All.GCaMP7s.animalID);
% gammaTable.CorrCoef = cat(1,data.CorrCoef.Rest.GCaMP7s.meanRs,data.CorrCoef.Whisk.GCaMP7s.meanRs,data.CorrCoef.NREM.GCaMP7s.meanRs,data.CorrCoef.REM.GCaMP7s.meanRs,...
%     data.CorrCoef.Awake.GCaMP7s.meanRs,data.CorrCoef.Sleep.GCaMP7s.meanRs,data.CorrCoef.All.GCaMP7s.meanRs);
% gammaTable.Behavior = cat(1,data.CorrCoef.Rest.GCaMP7s.behavior,data.CorrCoef.Whisk.GCaMP7s.behavior,data.CorrCoef.NREM.GCaMP7s.behavior,data.CorrCoef.REM.GCaMP7s.behavior,...
%     data.CorrCoef.Awake.GCaMP7s.behavior,data.CorrCoef.Sleep.GCaMP7s.behavior,data.CorrCoef.All.GCaMP7s.behavior);
% gammaFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
% gammaStats = fitglme(gammaTable,gammaFitFormula);
%% Fig. 7
summaryFigure = figure('Name','Fig7 (a-g)');
sgtitle('Figure 7 - Turner et al. 2020')
CC_xInds = ones(1,length(IOS_animalIDs));
CC_xInds2 = ones(1,length(data.CorrCoef.Awake.GCaMP7s.animalID));
% CC_xInds3 = ones(1,length(data.CorrCoef.Sleep.GCaMP7s.animalID));
%% [7a] power spectra of gamma-band power during different arousal-states
ax1 = subplot(3,3,1);
L1 = loglog(data.PowerSpec.Rest.GCaMP7s.meanCortf,data.PowerSpec.Rest.GCaMP7s.meanCortS,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.15,0.1 - 0.005,95],'FaceColor','w','EdgeColor','w')
L2 = loglog(data.PowerSpec.NREM.GCaMP7s.meanCortf,data.PowerSpec.NREM.GCaMP7s.meanCortS,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.15,1/30 - 0.005,95],'FaceColor','w','EdgeColor','w')
L3 = loglog(data.PowerSpec.REM.GCaMP7s.meanCortf,data.PowerSpec.REM.GCaMP7s.meanCortS,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.15,1/60 - 0.005,95],'FaceColor','w','EdgeColor','w')
L4 = loglog(data.PowerSpec.Awake.GCaMP7s.meanCortf,data.PowerSpec.Awake.GCaMP7s.meanCortS,'color',colorAlert,'LineWidth',2);
L5 = loglog(data.PowerSpec.Sleep.GCaMP7s.meanCortf,data.PowerSpec.Sleep.GCaMP7s.meanCortS,'color',colorAsleep,'LineWidth',2);
L6 = loglog(data.PowerSpec.All.GCaMP7s.meanCortf,data.PowerSpec.All.GCaMP7s.meanCortS,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title({'[7a] Cortical power','Gamma-band [30-100 Hz]'})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([L1,L2,L3,L4,L5,L6],'Rest','NREM','REM','Alert','Asleep','All','Location','SouthEast')
axis square
xlim([0.003,0.5])
ylim([0.1,100])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [7b] coherence^2 between bilateral gamma-band power during different arousal-states
ax2 = subplot(3,3,2);
semilogx(data.Coherr.Rest.GCaMP7s.meanf,data.Coherr.Rest.GCaMP7s.meanC.^2,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.05,0.1 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.Coherr.NREM.GCaMP7s.meanf,data.Coherr.NREM.GCaMP7s.meanC.^2,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/30 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.Coherr.REM.GCaMP7s.meanf,data.Coherr.REM.GCaMP7s.meanC.^2,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/60 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.Coherr.Awake.GCaMP7s.meanf,data.Coherr.Awake.GCaMP7s.meanC.^2,'color',colorAlert,'LineWidth',2);
% semilogx(data.Coherr.Sleep.GCaMP7s.meanf,data.Coherr.Sleep.GCaMP7s.meanC.^2,'color',colorAsleep,'LineWidth',2);
semilogx(data.Coherr.All.GCaMP7s.meanf,data.Coherr.All.GCaMP7s.meanC.^2,'color',colorAll,'LineWidth',2);
% confidence lines
% semilogx(data.Coherr.Rest.GCaMP7s.meanf,data.Coherr.Rest.GCaMP7s.maxConfC_Y.^2,'-','color',colorRest,'LineWidth',1);
% semilogx(data.Coherr.NREM.GCaMP7s.meanf,data.Coherr.NREM.GCaMP7s.maxConfC_Y.^2,'-','color',colorNREM,'LineWidth',1);
% semilogx(data.Coherr.REM.GCaMP7s.meanf,data.Coherr.REM.GCaMP7s.maxConfC_Y.^2,'-','color',colorREM,'LineWidth',1);
% semilogx(data.Coherr.Awake.GCaMP7s.meanf,data.Coherr.Awake.GCaMP7s.maxConfC_Y.^2,'-','color',colorAlert,'LineWidth',1);
% semilogx(data.Coherr.Sleep.GCaMP7s.meanf,data.Coherr.Sleep.GCaMP7s.maxConfC_Y.^2,'-','color',colorAsleep,'LineWidth',1);
% semilogx(data.Coherr.All.GCaMP7s.meanf,data.Coherr.All.GCaMP7s.maxConfC_Y.^2,'-','color',colorAll,'LineWidth',1);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'Bilateral coherence^2','GCaMP7s'})
axis square
xlim([0.01,0.5])
ylim([0,1])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [7c] Pearson's correlations between bilateral gamma-band power during different arousal-states
ax3 = subplot(3,3,3);
s1 = scatter(CC_xInds*1,data.CorrCoef.Rest.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.CorrCoef.Rest.GCaMP7s.meanR,data.CorrCoef.Rest.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(CC_xInds*2,data.CorrCoef.Whisk.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.CorrCoef.Whisk.GCaMP7s.meanR,data.CorrCoef.Whisk.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(CC_xInds*3,data.CorrCoef.NREM.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.CorrCoef.NREM.GCaMP7s.meanR,data.CorrCoef.NREM.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(CC_xInds*4,data.CorrCoef.REM.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.CorrCoef.REM.GCaMP7s.meanR,data.CorrCoef.REM.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
s5 = scatter(CC_xInds2*5,data.CorrCoef.Awake.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.CorrCoef.Awake.GCaMP7s.meanR,data.CorrCoef.Awake.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
% s6 = scatter(CC_xInds3*6,data.CorrCoef.Sleep.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
% e6 = errorbar(6,data.CorrCoef.Sleep.GCaMP7s.meanR,data.CorrCoef.Sleep.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e6.Color = 'black';
% e6.MarkerSize = 10;
% e6.CapSize = 10;
s7 = scatter(CC_xInds*7,data.CorrCoef.All.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e7 = errorbar(7,data.CorrCoef.All.GCaMP7s.meanR,data.CorrCoef.All.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
title({' Cortical Pearson''s corr. coef','GCaMP7s'})
ylabel('Corr. coefficient')
legend([s1,s2,s3,s4,s5,s7],'Rest','Whisk','NREM','REM','Alert','All','location','best')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields3) + 1])
% ylim([-0.1,1])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [7d] power spectra of GCaMP7s during different arousal-states
ax4 = subplot(3,3,4);
loglog(data.PowerSpec.Rest.GCaMP7s.meanCortf,data.PowerSpec.Rest.GCaMP7s.meanCortS,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.015,0.1 - 0.005,950],'FaceColor','w','EdgeColor','w')
loglog(data.PowerSpec.NREM.GCaMP7s.meanCortf,data.PowerSpec.NREM.GCaMP7s.meanCortS,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.015,1/30 - 0.005,950],'FaceColor','w','EdgeColor','w')
loglog(data.PowerSpec.REM.GCaMP7s.meanCortf,data.PowerSpec.REM.GCaMP7s.meanCortS,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.015,1/60 - 0.005,950],'FaceColor','w','EdgeColor','w')
loglog(data.PowerSpec.Awake.GCaMP7s.meanCortf,data.PowerSpec.Awake.GCaMP7s.meanCortS,'color',colorAlert,'LineWidth',2);
loglog(data.PowerSpec.Sleep.GCaMP7s.meanCortf,data.PowerSpec.Sleep.GCaMP7s.meanCortS,'color',colorAsleep,'LineWidth',2);
loglog(data.PowerSpec.All.GCaMP7s.meanCortf,data.PowerSpec.All.GCaMP7s.meanCortS,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title({'[7d] Cortical power','\Delta[GCaMP7s] (\muM)'})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.003,0.5]);
ylim([0.01,1000])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [7b] coherence^2 between bilateral GCaMP7s during different arousal-states
ax5 = subplot(3,3,5);
semilogx(data.Coherr.Rest.GCaMP7s.meanf,data.Coherr.Rest.GCaMP7s.meanC.^2,'color',colorRest,'LineWidth',2);
hold on
rectangle('Position',[0.005,0.05,0.1 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.Coherr.NREM.GCaMP7s.meanf,data.Coherr.NREM.GCaMP7s.meanC.^2,'color',colorNREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/30 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.Coherr.REM.GCaMP7s.meanf,data.Coherr.REM.GCaMP7s.meanC.^2,'color',colorREM,'LineWidth',2);
rectangle('Position',[0.005,0.05,1/60 - 0.005,0.90],'FaceColor','w','EdgeColor','w')
semilogx(data.Coherr.Awake.GCaMP7s.meanf,data.Coherr.Awake.GCaMP7s.meanC.^2,'color',colorAlert,'LineWidth',2);
% semilogx(data.Coherr.Sleep.GCaMP7s.meanf,data.Coherr.Sleep.GCaMP7s.meanC.^2,'color',colorAsleep,'LineWidth',2);
semilogx(data.Coherr.All.GCaMP7s.meanf,data.Coherr.All.GCaMP7s.meanC.^2,'color',colorAll,'LineWidth',2);
% confidence lines
% semilogx(data.Coherr.Rest.GCaMP7s.meanf,data.Coherr.Rest.GCaMP7s.maxConfC_Y.^2,'-','color',colorRest,'LineWidth',1);
% semilogx(data.Coherr.NREM.GCaMP7s.meanf,data.Coherr.NREM.GCaMP7s.maxConfC_Y.^2,'-','color',colorNREM,'LineWidth',1);
% semilogx(data.Coherr.REM.GCaMP7s.meanf,data.Coherr.REM.GCaMP7s.maxConfC_Y.^2,'-','color',colorREM,'LineWidth',1);
% semilogx(data.Coherr.Awake.GCaMP7s.meanf,data.Coherr.Awake.GCaMP7s.maxConfC_Y.^2,'-','color',colorAlert,'LineWidth',1);
% semilogx(data.Coherr.Sleep.GCaMP7s.meanf,data.Coherr.Sleep.GCaMP7s.maxConfC_Y.^2,'-','color',colorAsleep,'LineWidth',1);
% semilogx(data.Coherr.All.GCaMP7s.meanf,data.Coherr.All.GCaMP7s.maxConfC_Y.^2,'-','color',colorAll,'LineWidth',1);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'Bilateral coherence^2','\Delta[GCaMP7s]'})
axis square
xlim([0.01,0.5])
ylim([0,1])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [7e] Pearson's correlations between bilateral GCaMP7s during different arousal-states
ax6 = subplot(3,3,6);
scatter(CC_xInds*1,data.CorrCoef.Rest.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.CorrCoef.Rest.GCaMP7s.meanR,data.CorrCoef.Rest.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(CC_xInds*2,data.CorrCoef.Whisk.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.CorrCoef.Whisk.GCaMP7s.meanR,data.CorrCoef.Whisk.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(CC_xInds*3,data.CorrCoef.NREM.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.CorrCoef.NREM.GCaMP7s.meanR,data.CorrCoef.NREM.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(CC_xInds*4,data.CorrCoef.REM.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.CorrCoef.REM.GCaMP7s.meanR,data.CorrCoef.REM.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(CC_xInds2*5,data.CorrCoef.Awake.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAlert,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.CorrCoef.Awake.GCaMP7s.meanR,data.CorrCoef.Awake.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
% scatter(CC_xInds3*6,data.CorrCoef.Sleep.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAsleep,'jitter','on', 'jitterAmount',0.25);
% e6 = errorbar(6,data.CorrCoef.Sleep.GCaMP7s.meanR,data.CorrCoef.Sleep.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e6.Color = 'black';
% e6.MarkerSize = 10;
% e6.CapSize = 10;
scatter(CC_xInds*7,data.CorrCoef.All.GCaMP7s.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorAll,'jitter','on', 'jitterAmount',0.25);
e7 = errorbar(7,data.CorrCoef.All.GCaMP7s.meanR,data.CorrCoef.All.GCaMP7s.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
title({'Cortical Pearson''s corr. coef','\DeltaGCaMP7s'})
ylabel('Corr. coefficient')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields3) + 1])
% ylim([0,1])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];

%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig7']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig7'])
    %% statistical diary
%     diaryFile = [dirpath 'Fig7_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     % HbT statistical diary
%     disp('======================================================================================================================')
%     disp('[7c] Generalized linear mixed-effects model statistics for mean HbT corr. coef during Rest, Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(gammaStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Rest  Gamma P/P R: ' num2str(round(data.CorrCoef.Rest.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Rest.gammaBandPower.stdR,2))]); disp(' ')
%     disp(['Whisk Gamma P/P R: ' num2str(round(data.CorrCoef.Whisk.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Whisk.gammaBandPower.stdR,2))]); disp(' ')
%     disp(['NREM  Gamma P/P R: ' num2str(round(data.CorrCoef.NREM.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.NREM.gammaBandPower.stdR,2))]); disp(' ')
%     disp(['REM   Gamma P/P R: ' num2str(round(data.CorrCoef.REM.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.REM.gammaBandPower.stdR,2))]); disp(' ')
%     disp(['Awake Gamma P/P R: ' num2str(round(data.CorrCoef.Awake.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Awake.gammaBandPower.stdR,2))]); disp(' ')
%     disp(['Sleep Gamma P/P R: ' num2str(round(data.CorrCoef.Sleep.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Sleep.gammaBandPower.stdR,2))]); disp(' ')
%     disp(['All   Gamma P/P R: ' num2str(round(data.CorrCoef.All.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.All.gammaBandPower.stdR,2))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     % gamma statistical diary
%     disp('======================================================================================================================')
%     disp('[7f] Generalized linear mixed-effects model statistics for mean gamma-band corr. coef during Rest, Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(HbTStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Rest  [HbT] (uM) R: ' num2str(round(data.CorrCoef.Rest.HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Rest.HbT.stdR,2))]); disp(' ')
%     disp(['Whisk [HbT] (uM) R: ' num2str(round(data.CorrCoef.Whisk.HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Whisk.HbT.stdR,2))]); disp(' ')
%     disp(['NREM  [HbT] (uM) R: ' num2str(round(data.CorrCoef.NREM.HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.NREM.HbT.stdR,2))]); disp(' ')
%     disp(['REM   [HbT] (uM) R: ' num2str(round(data.CorrCoef.REM.HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.REM.HbT.stdR,2))]); disp(' ')
%     disp(['Awake [HbT] (uM) R: ' num2str(round(data.CorrCoef.Awake.HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Awake.HbT.stdR,2))]); disp(' ')
%     disp(['Sleep [HbT] (uM) R: ' num2str(round(data.CorrCoef.Sleep.HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Sleep.HbT.stdR,2))]); disp(' ')
%     disp(['All   [HbT] (uM) R: ' num2str(round(data.CorrCoef.All.HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.All.HbT.stdR,2))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     diary off
%     %% organized for supplemental table
%     % variable names
%     ColumnNames_R = {'Rest','Whisk','NREM','REM','Awake','All'};
%     % gamma-band R
%     for aa = 1:length(ColumnNames_R)
%         Gamma_R_MeanStD{1,aa} = [num2str(round(data.CorrCoef.(ColumnNames_R{1,aa}).gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.(ColumnNames_R{1,aa}).gammaBandPower.stdR,2))]; %#ok<*AGROW>
%     end
%     % gamma-band R p-values
% %     for aa = 1:length(ColumnNames_R)
% %         if strcmp(ColumnNames_R{1,aa},'Rest') == true
% %             Gamma_R_pVal{1,aa} = {' '};
% %         else
% %             Gamma_R_pVal{1,aa} = ['p < ' num2str(gammaStats.Coefficients.pValue(aa,1))];
% %         end
% %     end
%     % HbT R
%     for aa = 1:length(ColumnNames_R)
%         HbT_MeanStD{1,aa} = [num2str(round(data.CorrCoef.(ColumnNames_R{1,aa}).HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.(ColumnNames_R{1,aa}).HbT.stdR,2))];
%     end
%     % HbT R p-values
% %     for aa = 1:length(ColumnNames_R)
% %         if strcmp(ColumnNames_R{1,aa},'Rest') == true
% %             HbT_R_pVal{1,aa} = {' '};
% %         else
% %             HbT_R_pVal{1,aa} = ['p < ' num2str(HbTStats.Coefficients.pValue(aa,1))];
% %         end
% %     end
%     %% save table data
%     if isfield(AnalysisResults,'CorrCoef') == false
%         AnalysisResults.CorrCoef = [];
%     end
%     if isfield(AnalysisResults.CorrCoef,'gammaBandPower') == false
%         AnalysisResults.CorrCoef.columnNames = ColumnNames_R;
%         AnalysisResults.CorrCoef.gammaBandPower.meanStD = Gamma_R_MeanStD;
%         AnalysisResults.CorrCoef.gammaBandPower.p = Gamma_R_pVal;
%         AnalysisResults.CorrCoef.HbT.meanStD = HbT_MeanStD;
%         AnalysisResults.CorrCoef.HbT.p = HbT_R_pVal;
%         cd(rootFolder)
%         save('AnalysisResults.mat','AnalysisResults')
%     end
end

end
