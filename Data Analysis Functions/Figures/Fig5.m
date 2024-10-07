function [AnalysisResults] = Fig5(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 5 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
% colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
colorRest = [(0/256),(166/256),(81/256)];
colorWhisk = [(31/256),(120/256),(179/256)];
colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
IOS_animalIDs = {'T281','T282','T284','T285'};

%% mean TRITC comparison between behaviors
% pre-allocate the date for each day
IOS_behavFields = {'Rest','Whisk','Stim','NREM'};
for aa = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,aa};
    whiskFileIDs = unique(AnalysisResults.(animalID).MeanTRITC.Whisk.TRITC.FileIDs);
    whiskFileDates = [];
    % identify the unique days present for each animal using the whisking field.
    for bb = 1:length(whiskFileIDs)
        whiskFileDates{bb,1} = ConvertDate_IOS(whiskFileIDs{bb,1}); %#ok<*AGROW>
    end
    uniqueWhiskFileDates = unique(whiskFileDates);
    % put pre-allocate each date
    for dd = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,dd};
        for ee = 1:length(uniqueWhiskFileDates)
            fileDate = uniqueWhiskFileDates{ee,1};
            data.TRITC.(animalID).(behavField).(fileDate).MeanLH = [];
            data.TRITC.(animalID).(behavField).(fileDate).MeanRH = [];
            data.TRITC.(animalID).(behavField).(fileDate).IndLH = {};
            data.TRITC.(animalID).(behavField).(fileDate).IndRH = {};
            data.GCaMP7s.(animalID).(behavField).(fileDate).MeanLH = [];
            data.GCaMP7s.(animalID).(behavField).(fileDate).MeanRH = [];
            data.GCaMP7s.(animalID).(behavField).(fileDate).IndLH = {};
            data.GCaMP7s.(animalID).(behavField).(fileDate).IndRH = {};
        end
        procData.TRITC.(behavField).animalID{aa,1} = animalID;
        procData.TRITC.(behavField).behavior{aa,1} = behavField;
        procData.TRITC.(behavField).LH{aa,1} = 'LH';
        procData.TRITC.(behavField).RH{aa,1} = 'RH';
        procData.GCaMP7s.(behavField).animalID{aa,1} = animalID;
        procData.GCaMP7s.(behavField).behavior{aa,1} = behavField;
        procData.GCaMP7s.(behavField).LH{aa,1} = 'LH';
        procData.GCaMP7s.(behavField).RH{aa,1} = 'RH';
    end
end
%% TRITC
% put data into cell for each unique date
for ff = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,ff};
    for gg = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,gg};
        % data is structured slightly differently depending on class
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'Whisk') == true
            fileIDs = AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{hh,1});
                data.TRITC.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.TRITC.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.MeanLH(hh,1));
                data.TRITC.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.TRITC.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.MeanRH(hh,1));
                data.TRITC.(animalID).(behavField).(fileDate).IndLH = cat(1,data.TRITC.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndLH{hh,1});
                data.TRITC.(animalID).(behavField).(fileDate).IndRH = cat(1,data.TRITC.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndRH{hh,1});
            end
        elseif strcmp(behavField,'Stim') == true
            % left hem stims
            fileIDs = AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.LH_FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{hh,1});
                data.TRITC.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.TRITC.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.MeanLH(hh,1));
                data.TRITC.(animalID).(behavField).(fileDate).IndLH = cat(1,data.TRITC.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndLH{hh,1});
            end
            % right hem stims
            fileIDs = AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.RH_FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{hh,1});
                data.TRITC.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.TRITC.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.MeanRH(hh,1));
                data.TRITC.(animalID).(behavField).(fileDate).IndRH = cat(1,data.TRITC.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndRH{hh,1});
            end
        else
            fileIDs = AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.FileIDs;
            for ii = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{ii,1});
                data.TRITC.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.TRITC.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.MeanLH(ii,1));
                data.TRITC.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.TRITC.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.MeanRH(ii,1));
                data.TRITC.(animalID).(behavField).(fileDate).IndLH = cat(1,data.TRITC.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndLH{ii,1});
                data.TRITC.(animalID).(behavField).(fileDate).IndRH = cat(1,data.TRITC.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndRH{ii,1});
            end
        end
    end
end
% find the mean of the 10-second resting periods from each day to determine a baseline
for jj = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,jj};
    whiskFileIDs = unique(AnalysisResults.(animalID).MeanTRITC.Whisk.TRITC.FileIDs);
    whiskFileDates = [];
    for kk = 1:length(whiskFileIDs)
        whiskFileDates{kk,1} = ConvertDate_IOS(whiskFileIDs{kk,1}); %#ok<*AGROW>
    end
    uniqueWhiskFileDates = unique(whiskFileDates);
    % take mean from each day. Days with no data will show up as NaN and be excluded
    for ll = 1:length(uniqueWhiskFileDates)
        fileDate = uniqueWhiskFileDates{ll,1};
        data.TRITC.(animalID).Rest.(fileDate).baselineLH = mean(data.TRITC.(animalID).Rest.(fileDate).MeanLH);
        data.TRITC.(animalID).Rest.(fileDate).baselineRH = mean(data.TRITC.(animalID).Rest.(fileDate).MeanRH);
    end
end
% subtract the 10-second resting baseline for each day from the other data types. If the day doesn't have resting data,
% exclude it from analysis
for mm = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,mm};
    for nn = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,nn};
        % subtract each day's 10-second baseline from each behavior field
        fileDates = fieldnames(data.TRITC.(animalID).(behavField));
        for oo = 1:length(fileDates)
            fileDate = fileDates{oo,1};
            if strcmp(behavField,'Stim') == true || strcmp(behavField,'Whisk') == true
                data.TRITC.(animalID).(behavField).(fileDate).CorrMeanLH = data.TRITC.(animalID).(behavField).(fileDate).MeanLH;% - data.TRITC.(animalID).Rest.(fileDate).baselineLH;
                data.TRITC.(animalID).(behavField).(fileDate).CorrMeanRH = data.TRITC.(animalID).(behavField).(fileDate).MeanRH;%; - data.TRITC.(animalID).Rest.(fileDate).baselineRH;
                for pp = 1:length(data.TRITC.(animalID).(behavField).(fileDate).IndLH)
                    data.TRITC.(animalID).(behavField).(fileDate).CorrIndLH{pp,1} = data.TRITC.(animalID).(behavField).(fileDate).IndLH{pp,1};% - data.TRITC.(animalID).Rest.(fileDate).baselineLH;
                end
                for pp = 1:length(data.TRITC.(animalID).(behavField).(fileDate).IndRH)
                    data.TRITC.(animalID).(behavField).(fileDate).CorrIndRH{pp,1} = data.TRITC.(animalID).(behavField).(fileDate).IndRH{pp,1};% - data.TRITC.(animalID).Rest.(fileDate).baselineRH;
                end
            else
                data.TRITC.(animalID).(behavField).(fileDate).CorrMeanLH = data.TRITC.(animalID).(behavField).(fileDate).MeanLH - data.TRITC.(animalID).Rest.(fileDate).baselineLH;
                data.TRITC.(animalID).(behavField).(fileDate).CorrMeanRH = data.TRITC.(animalID).(behavField).(fileDate).MeanRH - data.TRITC.(animalID).Rest.(fileDate).baselineRH;
                for pp = 1:length(data.TRITC.(animalID).(behavField).(fileDate).IndLH)
                    data.TRITC.(animalID).(behavField).(fileDate).CorrIndLH{pp,1} = data.TRITC.(animalID).(behavField).(fileDate).IndLH{pp,1} - data.TRITC.(animalID).Rest.(fileDate).baselineLH;
                end
                for pp = 1:length(data.TRITC.(animalID).(behavField).(fileDate).IndRH)
                    data.TRITC.(animalID).(behavField).(fileDate).CorrIndRH{pp,1} = data.TRITC.(animalID).(behavField).(fileDate).IndRH{pp,1} - data.TRITC.(animalID).Rest.(fileDate).baselineRH;
                end
            end
        end
    end
end
% take the mean of the corrected data from each unique day
for qq = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,qq};
    for rr = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,rr};
        fileDates = fieldnames(data.TRITC.(animalID).(behavField));
        for ss = 1:length(fileDates)
            fileDate = fileDates{ss,1};
            data.TRITC.(animalID).(behavField).(fileDate).DayMeanLH = mean(data.TRITC.(animalID).(behavField).(fileDate).CorrMeanLH);
            data.TRITC.(animalID).(behavField).(fileDate).DayMeanRH = mean(data.TRITC.(animalID).(behavField).(fileDate).CorrMeanRH);
            data.TRITC.(animalID).(behavField).(fileDate).DayAllMeanLH = [];
            data.TRITC.(animalID).(behavField).(fileDate).DayAllMeanRH = [];
            data.TRITC.(animalID).(behavField).(fileDate).DayIndLH = [];
            data.TRITC.(animalID).(behavField).(fileDate).DayIndRH = [];
            % concatenate individual trials into a single array for each unique day
            if isfield(data.TRITC.(animalID).(behavField).(fileDate),'CorrIndLH') == true
                % left means - diff loop is necessary as STIM field has diff number of events
                for tt = 1:length(data.TRITC.(animalID).(behavField).(fileDate).CorrMeanLH)
                    data.TRITC.(animalID).(behavField).(fileDate).DayAllMeanLH = cat(2,data.TRITC.(animalID).(behavField).(fileDate).DayIndLH,data.TRITC.(animalID).(behavField).(fileDate).CorrIndLH{tt,1});
                end
                % right means
                for tt = 1:length(data.TRITC.(animalID).(behavField).(fileDate).CorrMeanRH)
                    data.TRITC.(animalID).(behavField).(fileDate).DayAllMeanRH = cat(2,data.TRITC.(animalID).(behavField).(fileDate).DayIndRH,data.TRITC.(animalID).(behavField).(fileDate).CorrIndRH{tt,1});
                end
                % left individual data pts
                for tt = 1:length(data.TRITC.(animalID).(behavField).(fileDate).CorrIndLH)
                    data.TRITC.(animalID).(behavField).(fileDate).DayIndLH = cat(2,data.TRITC.(animalID).(behavField).(fileDate).DayIndLH,data.TRITC.(animalID).(behavField).(fileDate).CorrIndLH{tt,1});
                end
                % right individual data pts
                for tt = 1:length(data.TRITC.(animalID).(behavField).(fileDate).CorrIndRH)
                    data.TRITC.(animalID).(behavField).(fileDate).DayIndRH = cat(2,data.TRITC.(animalID).(behavField).(fileDate).DayIndRH,data.TRITC.(animalID).(behavField).(fileDate).CorrIndRH{tt,1});
                end
            end
        end
    end
end
% put all the corrected means from each unique day into a single vector
nans = 1;
for uu = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,uu};
    for vv = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,vv};
        fileDates = fieldnames(data.TRITC.(animalID).(behavField));
        procData.TRITC.(animalID).(behavField).DayMeansLH = [];
        procData.TRITC.(animalID).(behavField).DayMeansRH = [];
        procData.TRITC.(animalID).(behavField).CatIndLH = [];
        procData.TRITC.(animalID).(behavField).CatIndRH = [];
        for ww = 1:length(fileDates)
            fileDate = fileDates{ww,1};
            if isnan(data.TRITC.(animalID).(behavField).(fileDate).DayMeanLH) == false
                procData.TRITC.(animalID).(behavField).DayMeansLH = cat(1,procData.TRITC.(animalID).(behavField).DayMeansLH,data.TRITC.(animalID).(behavField).(fileDate).DayMeanLH);
                procData.TRITC.(animalID).(behavField).DayMeansRH = cat(1,procData.TRITC.(animalID).(behavField).DayMeansRH,data.TRITC.(animalID).(behavField).(fileDate).DayMeanRH);
                procData.TRITC.(animalID).(behavField).CatIndLH = cat(2,procData.TRITC.(animalID).(behavField).CatIndLH,data.TRITC.(animalID).(behavField).(fileDate).DayIndLH);
                procData.TRITC.(animalID).(behavField).CatIndRH = cat(2,procData.TRITC.(animalID).(behavField).CatIndRH,data.TRITC.(animalID).(behavField).(fileDate).DayIndRH);
            else
                nans = nans + 1;
            end
        end
    end
end
% put all the means (of the corrected means) from each unique day into a single vector
for yy = 1:length(IOS_behavFields)
    behavField = IOS_behavFields{1,yy};
    procData.TRITC.(behavField).IndMeanTRITC = [];
    procData.TRITC.(behavField).CatTRITC = [];
    procData.TRITC.(behavField).meanLH = [];
    procData.TRITC.(behavField).meanRH = [];
    for zz = 1:length(IOS_animalIDs)
        animalID = IOS_animalIDs{1,zz};
        procData.TRITC.(behavField).IndMeanTRITC = cat(1,procData.TRITC.(behavField).IndMeanTRITC,mean(procData.TRITC.(animalID).(behavField).DayMeansLH),mean(procData.TRITC.(animalID).(behavField).DayMeansRH));
        procData.TRITC.(behavField).meanLH = cat(1,procData.TRITC.(behavField).meanLH,mean(procData.TRITC.(animalID).(behavField).DayMeansLH));
        procData.TRITC.(behavField).meanRH = cat(1,procData.TRITC.(behavField).meanRH,mean(procData.TRITC.(animalID).(behavField).DayMeansRH));
        procData.TRITC.(behavField).CatTRITC = cat(2,procData.TRITC.(behavField).CatTRITC,procData.TRITC.(animalID).(behavField).CatIndLH,procData.TRITC.(animalID).(behavField).CatIndRH);
    end
end
% take the mean and stdev across animals
for aaa = 1:length(IOS_behavFields)
    behavField = IOS_behavFields{1,aaa};
    procData.TRITC.(behavField).MeanTRITC = mean(procData.TRITC.(behavField).IndMeanTRITC,1);
    procData.TRITC.(behavField).StdMeanTRITC = std(procData.TRITC.(behavField).IndMeanTRITC,0,1);
end

%% GCaMP7s
        % put data into cell for each unique date
for ff = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,ff};
    for gg = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,gg};
        % data is structured slightly differently depending on class
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'Whisk') == true
            fileIDs = AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{hh,1});
                data.GCaMP7s.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.GCaMP7s.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.MeanLH(hh,1));
                data.GCaMP7s.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.GCaMP7s.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.MeanRH(hh,1));
                data.GCaMP7s.(animalID).(behavField).(fileDate).IndLH = cat(1,data.GCaMP7s.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndLH{hh,1});
                data.GCaMP7s.(animalID).(behavField).(fileDate).IndRH = cat(1,data.GCaMP7s.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndRH{hh,1});
            end
        elseif strcmp(behavField,'Stim') == true
            % left hem stims
            fileIDs = AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.LH_FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{hh,1});
                data.GCaMP7s.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.GCaMP7s.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.MeanLH(hh,1));
                data.GCaMP7s.(animalID).(behavField).(fileDate).IndLH = cat(1,data.GCaMP7s.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndLH{hh,1});
            end
            % right hem stims
            fileIDs = AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.RH_FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{hh,1});
                data.GCaMP7s.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.GCaMP7s.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.MeanRH(hh,1));
                data.GCaMP7s.(animalID).(behavField).(fileDate).IndRH = cat(1,data.GCaMP7s.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndRH{hh,1});
            end
        else
            fileIDs = AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.FileIDs;
            for ii = 1:length(fileIDs)
                fileDate = ConvertDate_IOS(fileIDs{ii,1});
                data.GCaMP7s.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.GCaMP7s.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.MeanLH(ii,1));
                data.GCaMP7s.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.GCaMP7s.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.MeanRH(ii,1));
                data.GCaMP7s.(animalID).(behavField).(fileDate).IndLH = cat(1,data.GCaMP7s.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndLH{ii,1});
                data.GCaMP7s.(animalID).(behavField).(fileDate).IndRH = cat(1,data.GCaMP7s.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndRH{ii,1});
            end
        end
    end
end
% find the mean of the 10-second resting periods from each day to determine a baseline
for jj = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,jj};
    whiskFileIDs = unique(AnalysisResults.(animalID).MeanGCaMP7s.Whisk.GCaMP7s.FileIDs);
    whiskFileDates = [];
    for kk = 1:length(whiskFileIDs)
        whiskFileDates{kk,1} = ConvertDate_IOS(whiskFileIDs{kk,1}); %#ok<*AGROW>
    end
    uniqueWhiskFileDates = unique(whiskFileDates);
    % take mean from each day. Days with no data will show up as NaN and be excluded
    for ll = 1:length(uniqueWhiskFileDates)
        fileDate = uniqueWhiskFileDates{ll,1};
        data.GCaMP7s.(animalID).Rest.(fileDate).baselineLH = mean(data.GCaMP7s.(animalID).Rest.(fileDate).MeanLH);
        data.GCaMP7s.(animalID).Rest.(fileDate).baselineRH = mean(data.GCaMP7s.(animalID).Rest.(fileDate).MeanRH);
    end
end
% subtract the 10-second resting baseline for each day from the other data types. If the day doesn't have resting data,
% exclude it from analysis
for mm = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,mm};
    for nn = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,nn};
        % subtract each day's 10-second baseline from each behavior field
        fileDates = fieldnames(data.GCaMP7s.(animalID).(behavField));
        for oo = 1:length(fileDates)
            fileDate = fileDates{oo,1};
            if strcmp(behavField,'Stim') == true || strcmp(behavField,'Whisk') == true
                data.GCaMP7s.(animalID).(behavField).(fileDate).CorrMeanLH = data.GCaMP7s.(animalID).(behavField).(fileDate).MeanLH;% - data.GCaMP7s.(animalID).Rest.(fileDate).baselineLH;
                data.GCaMP7s.(animalID).(behavField).(fileDate).CorrMeanRH = data.GCaMP7s.(animalID).(behavField).(fileDate).MeanRH;%; - data.GCaMP7s.(animalID).Rest.(fileDate).baselineRH;
                for pp = 1:length(data.GCaMP7s.(animalID).(behavField).(fileDate).IndLH)
                    data.GCaMP7s.(animalID).(behavField).(fileDate).CorrIndLH{pp,1} = data.GCaMP7s.(animalID).(behavField).(fileDate).IndLH{pp,1};% - data.GCaMP7s.(animalID).Rest.(fileDate).baselineLH;
                end
                for pp = 1:length(data.GCaMP7s.(animalID).(behavField).(fileDate).IndRH)
                    data.GCaMP7s.(animalID).(behavField).(fileDate).CorrIndRH{pp,1} = data.GCaMP7s.(animalID).(behavField).(fileDate).IndRH{pp,1};% - data.GCaMP7s.(animalID).Rest.(fileDate).baselineRH;
                end
            else
                data.GCaMP7s.(animalID).(behavField).(fileDate).CorrMeanLH = data.GCaMP7s.(animalID).(behavField).(fileDate).MeanLH - data.GCaMP7s.(animalID).Rest.(fileDate).baselineLH;
                data.GCaMP7s.(animalID).(behavField).(fileDate).CorrMeanRH = data.GCaMP7s.(animalID).(behavField).(fileDate).MeanRH - data.GCaMP7s.(animalID).Rest.(fileDate).baselineRH;
                for pp = 1:length(data.GCaMP7s.(animalID).(behavField).(fileDate).IndLH)
                    data.GCaMP7s.(animalID).(behavField).(fileDate).CorrIndLH{pp,1} = data.GCaMP7s.(animalID).(behavField).(fileDate).IndLH{pp,1} - data.GCaMP7s.(animalID).Rest.(fileDate).baselineLH;
                end
                for pp = 1:length(data.GCaMP7s.(animalID).(behavField).(fileDate).IndRH)
                    data.GCaMP7s.(animalID).(behavField).(fileDate).CorrIndRH{pp,1} = data.GCaMP7s.(animalID).(behavField).(fileDate).IndRH{pp,1} - data.GCaMP7s.(animalID).Rest.(fileDate).baselineRH;
                end
            end
        end
    end
end
% take the mean of the corrected data from each unique day
for qq = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,qq};
    for rr = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,rr};
        fileDates = fieldnames(data.GCaMP7s.(animalID).(behavField));
        for ss = 1:length(fileDates)
            fileDate = fileDates{ss,1};
            data.GCaMP7s.(animalID).(behavField).(fileDate).DayMeanLH = mean(data.GCaMP7s.(animalID).(behavField).(fileDate).CorrMeanLH);
            data.GCaMP7s.(animalID).(behavField).(fileDate).DayMeanRH = mean(data.GCaMP7s.(animalID).(behavField).(fileDate).CorrMeanRH);
            data.GCaMP7s.(animalID).(behavField).(fileDate).DayAllMeanLH = [];
            data.GCaMP7s.(animalID).(behavField).(fileDate).DayAllMeanRH = [];
            data.GCaMP7s.(animalID).(behavField).(fileDate).DayIndLH = [];
            data.GCaMP7s.(animalID).(behavField).(fileDate).DayIndRH = [];
            % concatenate individual trials into a single array for each unique day
            if isfield(data.GCaMP7s.(animalID).(behavField).(fileDate),'CorrIndLH') == true
                % left means - diff loop is necessary as STIM field has diff number of events
                for tt = 1:length(data.GCaMP7s.(animalID).(behavField).(fileDate).CorrMeanLH)
                    data.GCaMP7s.(animalID).(behavField).(fileDate).DayAllMeanLH = cat(2,data.GCaMP7s.(animalID).(behavField).(fileDate).DayIndLH,data.GCaMP7s.(animalID).(behavField).(fileDate).CorrIndLH{tt,1});
                end
                % right means
                for tt = 1:length(data.GCaMP7s.(animalID).(behavField).(fileDate).CorrMeanRH)
                    data.GCaMP7s.(animalID).(behavField).(fileDate).DayAllMeanRH = cat(2,data.GCaMP7s.(animalID).(behavField).(fileDate).DayIndRH,data.GCaMP7s.(animalID).(behavField).(fileDate).CorrIndRH{tt,1});
                end
                % left individual data pts
                for tt = 1:length(data.GCaMP7s.(animalID).(behavField).(fileDate).CorrIndLH)
                    data.GCaMP7s.(animalID).(behavField).(fileDate).DayIndLH = cat(2,data.GCaMP7s.(animalID).(behavField).(fileDate).DayIndLH,data.GCaMP7s.(animalID).(behavField).(fileDate).CorrIndLH{tt,1});
                end
                % right individual data pts
                for tt = 1:length(data.GCaMP7s.(animalID).(behavField).(fileDate).CorrIndRH)
                    data.GCaMP7s.(animalID).(behavField).(fileDate).DayIndRH = cat(2,data.GCaMP7s.(animalID).(behavField).(fileDate).DayIndRH,data.GCaMP7s.(animalID).(behavField).(fileDate).CorrIndRH{tt,1});
                end
            end
        end
    end
end
% put all the corrected means from each unique day into a single vector
nans = 1;
for uu = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,uu};
    for vv = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,vv};
        fileDates = fieldnames(data.GCaMP7s.(animalID).(behavField));
        procData.GCaMP7s.(animalID).(behavField).DayMeansLH = [];
        procData.GCaMP7s.(animalID).(behavField).DayMeansRH = [];
        procData.GCaMP7s.(animalID).(behavField).CatIndLH = [];
        procData.GCaMP7s.(animalID).(behavField).CatIndRH = [];
        for ww = 1:length(fileDates)
            fileDate = fileDates{ww,1};
            if isnan(data.GCaMP7s.(animalID).(behavField).(fileDate).DayMeanLH) == false
                procData.GCaMP7s.(animalID).(behavField).DayMeansLH = cat(1,procData.GCaMP7s.(animalID).(behavField).DayMeansLH,data.GCaMP7s.(animalID).(behavField).(fileDate).DayMeanLH);
                procData.GCaMP7s.(animalID).(behavField).DayMeansRH = cat(1,procData.GCaMP7s.(animalID).(behavField).DayMeansRH,data.GCaMP7s.(animalID).(behavField).(fileDate).DayMeanRH);
                procData.GCaMP7s.(animalID).(behavField).CatIndLH = cat(2,procData.GCaMP7s.(animalID).(behavField).CatIndLH,data.GCaMP7s.(animalID).(behavField).(fileDate).DayIndLH);
                procData.GCaMP7s.(animalID).(behavField).CatIndRH = cat(2,procData.GCaMP7s.(animalID).(behavField).CatIndRH,data.GCaMP7s.(animalID).(behavField).(fileDate).DayIndRH);
            else
                nans = nans + 1;
            end
        end
    end
end
% put all the means (of the corrected means) from each unique day into a single vector
for yy = 1:length(IOS_behavFields)
    behavField = IOS_behavFields{1,yy};
    procData.GCaMP7s.(behavField).IndMeanGCaMP7s = [];
    procData.GCaMP7s.(behavField).CatGCaMP7s = [];
    procData.GCaMP7s.(behavField).meanLH = [];
    procData.GCaMP7s.(behavField).meanRH = [];
    for zz = 1:length(IOS_animalIDs)
        animalID = IOS_animalIDs{1,zz};
        procData.GCaMP7s.(behavField).IndMeanGCaMP7s = cat(1,procData.GCaMP7s.(behavField).IndMeanGCaMP7s,mean(procData.GCaMP7s.(animalID).(behavField).DayMeansLH),mean(procData.GCaMP7s.(animalID).(behavField).DayMeansRH));
        procData.GCaMP7s.(behavField).meanLH = cat(1,procData.GCaMP7s.(behavField).meanLH,mean(procData.GCaMP7s.(animalID).(behavField).DayMeansLH));
        procData.GCaMP7s.(behavField).meanRH = cat(1,procData.GCaMP7s.(behavField).meanRH,mean(procData.GCaMP7s.(animalID).(behavField).DayMeansRH));
        procData.GCaMP7s.(behavField).CatGCaMP7s = cat(2,procData.GCaMP7s.(behavField).CatGCaMP7s,procData.GCaMP7s.(animalID).(behavField).CatIndLH,procData.GCaMP7s.(animalID).(behavField).CatIndRH);
    end
end
% take the mean and stdev across animals
for aaa = 1:length(IOS_behavFields)
    behavField = IOS_behavFields{1,aaa};
    procData.GCaMP7s.(behavField).MeanGCaMP7s = mean(procData.GCaMP7s.(behavField).IndMeanGCaMP7s,1);
    procData.GCaMP7s.(behavField).StdMeanGCaMP7s = std(procData.GCaMP7s.(behavField).IndMeanGCaMP7s,0,1);
end
%% statistics - generalized linear mixed effects model
% TRITCtableSize = cat(1,procData.TRITC.Rest.meanLH,procData.TRITC.Rest.meanRH,procData.TRITC.Whisk.meanLH,procData.TRITC.Whisk.meanRH,...
%     procData.TRITC.Stim.meanLH,procData.TRITC.Stim.meanRH,procData.TRITC.NREM.meanLH,procData.TRITC.NREM.meanRH,procData.TRITC.REM.meanLH,procData.TRITC.REM.meanRH);
% TRITCTable = table('Size',[size(TRITCtableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','TRITC','Behavior','Hemisphere'});
% TRITCTable.Mouse = cat(1,procData.TRITC.Rest.animalID,procData.TRITC.Rest.animalID,procData.TRITC.Whisk.animalID,procData.TRITC.Whisk.animalID,...
%     procData.TRITC.Stim.animalID,procData.TRITC.Stim.animalID,procData.TRITC.NREM.animalID,procData.TRITC.NREM.animalID,procData.TRITC.REM.animalID,procData.TRITC.REM.animalID);
% TRITCTable.TRITC = cat(1,procData.TRITC.Rest.meanLH,procData.TRITC.Rest.meanRH,procData.TRITC.Whisk.meanLH,procData.TRITC.Whisk.meanRH,...
%     procData.TRITC.Stim.meanLH,procData.TRITC.Stim.meanRH,procData.TRITC.NREM.meanLH,procData.TRITC.NREM.meanRH,procData.TRITC.REM.meanLH,procData.TRITC.REM.meanRH);
% TRITCTable.Behavior = cat(1,procData.TRITC.Rest.behavior,procData.TRITC.Rest.behavior,procData.TRITC.Whisk.behavior,procData.TRITC.Whisk.behavior,...
%     procData.TRITC.Stim.behavior,procData.TRITC.Stim.behavior,procData.TRITC.NREM.behavior,procData.TRITC.NREM.behavior,procData.TRITC.REM.behavior,procData.TRITC.REM.behavior);
% TRITCTable.Hemisphere = cat(1,procData.TRITC.Rest.LH,procData.TRITC.Rest.RH,procData.TRITC.Whisk.LH,procData.TRITC.Whisk.RH,...
%     procData.TRITC.Stim.LH,procData.TRITC.Stim.RH,procData.TRITC.NREM.LH,procData.TRITC.NREM.RH,procData.TRITC.REM.LH,procData.TRITC.REM.RH);
% TRITCFitFormula = 'TRITC ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
% TRITCStats = fitglme(TRITCTable,TRITCFitFormula); %#ok<*NASGU>

%% Fig. 5
summaryFigure = figure('Name','Fig5 (a-f)');
sgtitle('Figure 5 - Turner et al. 2020')
%% [5a] mean TRITC during different behaviors
ax1 = subplot(2,2,1);
TRITC_xInds = ones(1,length(IOS_animalIDs)*2);
s1 = scatter(TRITC_xInds*1,procData.TRITC.Rest.IndMeanTRITC,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,procData.TRITC.Rest.MeanTRITC,procData.TRITC.Rest.StdMeanTRITC,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(TRITC_xInds*2,procData.TRITC.Whisk.IndMeanTRITC,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,procData.TRITC.Whisk.MeanTRITC,procData.TRITC.Whisk.StdMeanTRITC,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(TRITC_xInds*3,procData.TRITC.Stim.IndMeanTRITC,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,procData.TRITC.Stim.MeanTRITC,procData.TRITC.Stim.StdMeanTRITC,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(TRITC_xInds*4,procData.TRITC.NREM.IndMeanTRITC,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,procData.TRITC.NREM.MeanTRITC,procData.TRITC.NREM.StdMeanTRITC,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
% s5 = scatter(TRITC_xInds*5,procData.TRITC.REM.IndMeanTRITC,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
% e5 = errorbar(5,procData.TRITC.REM.MeanTRITC,procData.TRITC.REM.StdMeanTRITC,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e5.Color = 'black';
% e5.MarkerSize = 10;
% e5.CapSize = 10;
title({'Mean \Delta[Rhodamine] ','during arousal-states'})
ylabel('\DeltaTRITC ')
legend([s1,s2,s3,s4],'Rest','Whisk','Stim','NREM','Location','best')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(IOS_behavFields) + 1])
% ylim([-10,100])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [5a bottom] mean TRITC distribution during different behaviors
ax3 = subplot(2,2,3);
[xCurve1,yCurve1] = SmoothHistogramBinsFit(procData.TRITC.Rest.CatTRITC,11,'normal');
[xCurve2,yCurve2] = SmoothHistogramBinsFit(procData.TRITC.Whisk.CatTRITC,11,'normal');
[xCurve3,yCurve3] = SmoothHistogramBinsFit(procData.TRITC.Stim.CatTRITC,11,'normal');
[xCurve4,yCurve4] = SmoothHistogramBinsFit(procData.TRITC.NREM.CatTRITC,11,'normal');
% [xCurve5,yCurve5] = SmoothHistogramBinsFit(procData.TRITC.REM.CatTRITC,11,'normal');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
% plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
title({'\Delta[Rhodamine] ','arousal-state distribution'})
xlabel('\Delta[Rhodamine] ')
ylabel('Probability')
axis square
set(gca,'box','off')
% xlim([-50,150])
% ylim([0,0.6])
ax3.TickLength = [0.03,0.03];
%% [5a] mean GCaMP7s during different behaviors
ax2 = subplot(2,2,2);
GCaMP7s_xInds = ones(1,length(IOS_animalIDs)*2);
s1 = scatter(GCaMP7s_xInds*1,procData.GCaMP7s.Rest.IndMeanGCaMP7s,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,procData.GCaMP7s.Rest.MeanGCaMP7s,procData.GCaMP7s.Rest.StdMeanGCaMP7s,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(GCaMP7s_xInds*2,procData.GCaMP7s.Whisk.IndMeanGCaMP7s,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,procData.GCaMP7s.Whisk.MeanGCaMP7s,procData.GCaMP7s.Whisk.StdMeanGCaMP7s,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(GCaMP7s_xInds*3,procData.GCaMP7s.Stim.IndMeanGCaMP7s,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,procData.GCaMP7s.Stim.MeanGCaMP7s,procData.GCaMP7s.Stim.StdMeanGCaMP7s,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(GCaMP7s_xInds*4,procData.GCaMP7s.NREM.IndMeanGCaMP7s,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,procData.GCaMP7s.NREM.MeanGCaMP7s,procData.GCaMP7s.NREM.StdMeanGCaMP7s,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
% s5 = scatter(GCaMP7s_xInds*5,procData.GCaMP7s.REM.IndMeanGCaMP7s,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
% e5 = errorbar(5,procData.GCaMP7s.REM.MeanGCaMP7s,procData.GCaMP7s.REM.StdMeanGCaMP7s,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e5.Color = 'black';
% e5.MarkerSize = 10;
% e5.CapSize = 10;
title({'Mean \Delta[GCaMP7s] ','during arousal-states'})
ylabel('\DeltaGCaMP7s ')
% legend([s1,s2,s3,s4,s5],'Rest','Whisk','Stim','NREM','REM','Location','bestoutside')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(IOS_behavFields) + 1])
% ylim([-10,100])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [5a bottom] mean GCaMP7s distribution during different behaviors
ax4 = subplot(2,2,4);
[xCurve1,yCurve1] = SmoothHistogramBinsFit(procData.GCaMP7s.Rest.CatGCaMP7s,11,'normal');
[xCurve2,yCurve2] = SmoothHistogramBinsFit(procData.GCaMP7s.Whisk.CatGCaMP7s,11,'normal');
[xCurve3,yCurve3] = SmoothHistogramBinsFit(procData.GCaMP7s.Stim.CatGCaMP7s,11,'normal');
[xCurve4,yCurve4] = SmoothHistogramBinsFit(procData.GCaMP7s.NREM.CatGCaMP7s,11,'normal');
% [xCurve5,yCurve5] = SmoothHistogramBinsFit(procData.GCaMP7s.REM.CatGCaMP7s,11,'normal');
plot(xCurve1,yCurve1,'color',colorRest,'LineWidth',2)
hold on
plot(xCurve2,yCurve2,'color',colorWhisk,'LineWidth',2)
plot(xCurve3,yCurve3,'color',colorStim,'LineWidth',2)
plot(xCurve4,yCurve4,'color',colorNREM,'LineWidth',2)
% plot(xCurve5,yCurve5,'color',colorREM,'LineWidth',2)
title({'\Delta[GCaMP7s] ','arousal-state distribution'})
xlabel('\Delta[GCaMP7s] ')
ylabel('Probability')
axis square
set(gca,'box','off')
% xlim([-50,150])
% ylim([0,0.6])
ax4.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    set(summaryFigure,'PaperPositionMode','auto');
    savefig(summaryFigure,[dirpath 'Fig5']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig5'])
%     %% statistical diary
%     diaryFile = [dirpath 'Fig5_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     % TRITC statistical diary
%     disp('======================================================================================================================')
%     disp('[5a] Generalized linear mixed-effects model statistics for mean TRITC during Rest, Whisk, Stim, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(TRITCStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Rest  [TRITC] (uM): Set = 0']); disp(' ') %#ok<NBRAK>
%     disp(['Whisk [TRITC] (uM): ' num2str(round(procData.TRITC.Whisk.MeanTRITC,1)) ' +/- ' num2str(round(procData.TRITC.Whisk.StdMeanTRITC,1))]); disp(' ')
%     disp(['Stim  [TRITC] (uM): ' num2str(round(procData.TRITC.Stim.MeanTRITC,1)) ' +/- ' num2str(round(procData.TRITC.Stim.StdMeanTRITC,1))]); disp(' ')
%     disp(['NREM  [TRITC] (uM): ' num2str(round(procData.TRITC.NREM.MeanTRITC,1)) ' +/- ' num2str(round(procData.TRITC.NREM.StdMeanTRITC,1))]); disp(' ')
%     disp(['REM   [TRITC] (uM): ' num2str(round(procData.TRITC.REM.MeanTRITC,1)) ' +/- ' num2str(round(procData.TRITC.REM.StdMeanTRITC,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     % Peak vessel diameter statistical diary
%     disp('======================================================================================================================')
%     disp('[5b] Generalized linear mixed-effects model statistics for mean vessel diameter during Rest, Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(vesselStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Rest  D/D (%): Set = 0']); disp(' ') %#ok<NBRAK>
%     disp(['Whisk D/D (%): ' num2str(round(procData.TwoP.Whisk.MeanDiam,1)) ' +/- ' num2str(round(procData.TwoP.Whisk.StdMeanDiam,1))]); disp(' ')
%     disp(['NREM  D/D (%): ' num2str(round(procData.TwoP.NREM.MeanDiam,1)) ' +/- ' num2str(round(procData.TwoP.NREM.StdMeanDiam,1))]); disp(' ')
%     disp(['REM   D/D (%): ' num2str(round(procData.TwoP.REM.MeanDiam,1)) ' +/- ' num2str(round(procData.TwoP.REM.StdMeanDiam,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     % LDF flow statistical diary
%     disp('======================================================================================================================')
%     disp('[5c] Generalized linear mixed-effects model statistics for mean doppler flow during Rest, Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(flowStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Rest  Q/Q (%): Set = 0']); disp(' ') %#ok<NBRAK>
%     disp(['Whisk Q/Q (%): ' num2str(round(procData.LDF.Whisk.MeanLDF,1)) ' +/- ' num2str(round(procData.LDF.Whisk.StdLDF,1))]); disp(' ')
%     disp(['NREM  Q/Q (%): ' num2str(round(procData.LDF.NREM.MeanLDF,1)) ' +/- ' num2str(round(procData.LDF.NREM.StdLDF,1))]); disp(' ')
%     disp(['REM   Q/Q (%): ' num2str(round(procData.LDF.REM.MeanLDF,1)) ' +/- ' num2str(round(procData.LDF.REM.StdLDF,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     diary off
end

end
