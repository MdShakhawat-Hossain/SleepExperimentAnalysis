  








%% apply TDT correction to the data and trim excess time
    DataTrim = floor(trimTime*tdtSamplingRate);
    tdtFields = {'RH','LH'};
    for cc = 1:length(tdtFields)
        subfields = fieldnames(FiberData.(tdtFields{1,cc}))';
        for dd = 1:length(subfields)
            subsubfields = fieldnames(FiberData.(tdtFields{1,cc}).(subfields{1,dd}))';
            for nn = 1:length( subsubfields)
                RawData.data.(char(tdtFields(1,cc))).(char(subfields(1,dd))).(char(subsubfields(1,nn))) = FiberData.(char(tdtFields(1,cc))).(char(subfields(1,dd))).(char(subsubfields(1,nn)))(floor(trimTime*tdtSamplingRate):end - (DataTrim + 1));
            end
        end
    end

    labviewFields = {'cortical_LH','cortical_RH','hippocampus','EMG','forceSensor','stimulations'};
    for cc = 1:length(labviewFields)
        labviewAnalogShift = horzcat(labviewAnalogPad,RawData.data.([labviewFields{1,cc} '_long']));
        RawData.data.(labviewFields{1,cc}) = labviewAnalogShift(floor(trimTime*analogSamplingRate):end - (labviewAnalogCut + 1));
    end