function [DataStruct] = NormBehavioralDataStruct_FP_EEG(DataStruct,RestingBaselines,baselineType)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%
% Originally written by Aaron T. Winder
%________________________________________________________________________________________________________________________
%
%   Purpose: Normalizes data structures based on resting values
%________________________________________________________________________________________________________________________

dataTypes = fieldnames(DataStruct);
for dT = 1:length(dataTypes)
    dataType = char(dataTypes(dT));
    if strcmp(dataType,'CBV') == false && strcmp(dataType,'GFP') == false && strcmp(dataType,'Isos') == false % normalize data for everything else
        hemisphereDataTypes = fieldnames(DataStruct.(dataType));       
        for hDT = 1:length(hemisphereDataTypes)
            hemDataType = char(hemisphereDataTypes(hDT));           
            if isfield(DataStruct.(dataType).(hemDataType),'whisk')
                behaviorFields = fieldnames(DataStruct.(dataType).(hemDataType));       
                for bF = 1:length(behaviorFields)
                    behField = char(behaviorFields(bF));
                    if ~isempty(DataStruct.(dataType).(hemDataType).(behField).data)
                        NormData = DataStruct.(dataType).(hemDataType).(behField).data;
                        [uniqueDays,~,~] = GetUniqueDays_FP(DataStruct.(dataType).(hemDataType).(behField).fileDates);
                        for uD = 1:length(uniqueDays)
                            date = uniqueDays{uD};
                            strDay = ConvertDate_FP(date);
                            [~,dayInds] = GetDayInds_FP(DataStruct.(dataType).(hemDataType).(behField).fileDates,date);                          
                            disp(['Normalizing ' (hemDataType) ' ' (dataType) ' ' (behField) ' for ' (strDay) '...']); disp(' ')
                            % Calculate the baseline differently depending on data type
                            if iscell(DataStruct.(dataType).(hemDataType).(behField).data)
                                dayData = DataStruct.(dataType).(hemDataType).(behField).data(dayInds);
                                normDayData = cell(size(dayData));
                                dayBaselineMean = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).mean;                               
                                for dD = 1:size(dayData,1)
                                    cellBaseMean = dayBaselineMean*ones(1,size(dayData{dD},2));
                                    normDayData{dD} = dayData{dD}./cellBaseMean - 1;
                                end
                                NormData(dayInds) = normDayData;
                            else
                                dayBaselineMean = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).mean;
                                % Preallocate array and use for permutation
                                normDayData = DataStruct.(dataType).(hemDataType).(behField).data(dayInds,:,:);
                                % Permute norm_session_data to handle both matrix and array (squeeze
                                % causes a matrix dimension error if not permuted)
                                dayData = permute(normDayData,unique([2,1,ndims(normDayData)],'stable'));
                                for dD = 1:size(dayData,2)
                                    normDayData(dD,:,:) = squeeze(dayData(:,dD,:))./(ones(size(dayData,1),1)*dayBaselineMean) - 1;
                                end
                                NormData(dayInds,:,:) = normDayData;
                            end
                            DataStruct.(dataType).(hemDataType).(behField).NormData = NormData;
                        end
                    end
                end
            else
                NormData = DataStruct.(dataType).(hemDataType).data;
                [uniqueDays,~,~] = GetUniqueDays_FP(DataStruct.(dataType).(hemDataType).fileDates);
                for uD = 1:length(uniqueDays)
                    date = uniqueDays{uD};
                    strDay = ConvertDate_FP(date);
                    [~,dayInds] = GetDayInds_FP(DataStruct.(dataType).(hemDataType).fileDates,date);             
                    disp(['Normalizing ' (hemDataType) ' ' (dataType) ' for ' (strDay) '...']); disp(' ')
                    % Calculate the baseline differently depending on data type
                    if iscell(DataStruct.(dataType).(hemDataType).data)
                        dayData = DataStruct.(dataType).(hemDataType).data(dayInds);
                        normDayData = cell(size(dayData));
                        dayBaselineMean = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).mean;
                        for dD = 1:size(dayData,1)
                            cellBaseMean = dayBaselineMean*ones(1,size(dayData{dD},2));
                            normDayData{dD} = dayData{dD}./cellBaseMean - 1;
                        end
                        NormData(dayInds) = normDayData;
                    else
                        dayBaselineMean = RestingBaselines.(dataType).(hemDataType).(strDay).mean;
                        % Preallocate array and use for permutation
                        normDayData = DataStruct.(dataType).(hemDataType).data(dayInds,:,:);
                        % Permute norm_session_data to handle both matrix and array (squeeze
                        % causes a matrix dimension error if not permuted)
                        dayData = permute(normDayData,unique([2,1,ndims(normDayData)],'stable'));
                        
                        for dD = 1:size(dayData,2)
                            normDayData(dD,:,:) = squeeze(dayData(:,dD,:))./(ones(size(dayData,1),1)*dayBaselineMean) - 1;
                        end
                        NormData(dayInds,:,:) = normDayData;
                    end
                    DataStruct.(dataType).(hemDataType).NormData = NormData;
                end
            end
        end
    elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true %% subtract the resting value the data for fiber photometry
        hemisphereDataTypes = fieldnames(DataStruct.(dataType));       
        for hDT = 1:length(hemisphereDataTypes)
            hemDataType = char(hemisphereDataTypes(hDT));           
            if isfield(DataStruct.(dataType).(hemDataType),'whisk')
                behaviorFields = fieldnames(DataStruct.(dataType).(hemDataType));       
                for bF = 1:length(behaviorFields)
                    behField = char(behaviorFields(bF));
                    if ~isempty(DataStruct.(dataType).(hemDataType).(behField).data)
                        NormData = DataStruct.(dataType).(hemDataType).(behField).data;
                        [uniqueDays,~,~] = GetUniqueDays_FP(DataStruct.(dataType).(hemDataType).(behField).fileDates);
                        for uD = 1:length(uniqueDays)
                            date = uniqueDays{uD};
                            strDay = ConvertDate_FP(date);
                            [~,dayInds] = GetDayInds_FP(DataStruct.(dataType).(hemDataType).(behField).fileDates,date);                          
                            disp(['Normalizing ' (hemDataType) ' ' (dataType) ' ' (behField) ' for ' (strDay) '...']); disp(' ')
                            % Calculate the baseline differently depending on data type
                            if iscell(DataStruct.(dataType).(hemDataType).(behField).data)
                                dayData = DataStruct.(dataType).(hemDataType).(behField).data(dayInds);
                                normDayData = cell(size(dayData));
                                dayBaselineMean = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).mean;
                                dayBaselineStd = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).std;                               

                                for dD = 1:size(dayData,1)
                                    cellBaseMean = dayBaselineMean*ones(1,size(dayData{dD},2));
                                    cellBaseStd = dayBaselineStd*ones(1,size(dayData{dD},2));
                                    normDayData{dD} = (dayData{dD}-cellBaseMean)./cellBaseStd;
                                end
                                NormData(dayInds) = normDayData;
                            else
                                dayBaselineMean = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).mean;
                                dayBaselineStd = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).std;                               
                                % Preallocate array and use for permutation
                                normDayData = DataStruct.(dataType).(hemDataType).(behField).data(dayInds,:,:);
                                % Permute norm_session_data to handle both matrix and array (squeeze
                                % causes a matrix dimension error if not permuted)
                                dayData = permute(normDayData,unique([2,1,ndims(normDayData)],'stable'));
                                for dD = 1:size(dayData,2)
                                    normDayData(dD,:,:) = (squeeze(dayData(:,dD,:))-(ones(size(dayData,1),1)*dayBaselineMean));%./ (ones(size(dayData,1),1)*dayBaselineStd);
                                end
                                NormData(dayInds,:,:) = normDayData;
                            end
                            DataStruct.(dataType).(hemDataType).(behField).NormData = NormData;
                        end
                    end
                end
            else
                NormData = DataStruct.(dataType).(hemDataType).data;
                [uniqueDays,~,~] = GetUniqueDays_FP(DataStruct.(dataType).(hemDataType).fileDates);
                for uD = 1:length(uniqueDays)
                    date = uniqueDays{uD};
                    strDay = ConvertDate_FP(date);
                    [~,dayInds] = GetDayInds_FP(DataStruct.(dataType).(hemDataType).fileDates,date);             
                    disp(['Normalizing ' (hemDataType) ' ' (dataType) ' for ' (strDay) '...']); disp(' ')
                    % Calculate the baseline differently depending on data type
                    if iscell(DataStruct.(dataType).(hemDataType).data)
                        dayData = DataStruct.(dataType).(hemDataType).data(dayInds);
                        normDayData = cell(size(dayData));
                        dayBaselineMean = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).mean;
                        dayBaselineStd = RestingBaselines.(baselineType).(dataType).(hemDataType).(strDay).std;
                        for dD = 1:size(dayData,1)
                            cellBaseMean = dayBaselineMean*ones(1,size(dayData{dD},2));
                            cellBaseStd = dayBaselineStd*ones(1,size(dayData{dD},2));
                            normDayData{dD} = (dayData{dD}-cellBaseMean)./cellBaseStd;
                        end
                        NormData(dayInds) = normDayData;
                    else
                        dayBaselineMean = RestingBaselines.(dataType).(hemDataType).(strDay).mean;
                        dayBaselineStd = RestingBaselines.(dataType).(hemDataType).(strDay).std;
                        % Preallocate array and use for permutation
                        normDayData = DataStruct.(dataType).(hemDataType).data(dayInds,:,:);
                        % Permute norm_session_data to handle both matrix and array (squeeze
                        % causes a matrix dimension error if not permuted)
                        dayData = permute(normDayData,unique([2,1,ndims(normDayData)],'stable'));
                        
                        for dD = 1:size(dayData,2)
                            normDayData(dD,:,:) = (squeeze(dayData(:,dD,:))-(ones(size(dayData,1),1)*dayBaselineMean));%./ (ones(size(dayData,1),1)*dayBaselineStd);
                        end
                        NormData(dayInds,:,:) = normDayData;
                    end
                    DataStruct.(dataType).(hemDataType).NormData = NormData;
                end
            end
        end
    end

end
