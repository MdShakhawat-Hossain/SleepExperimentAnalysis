function UpdateMicroArousals(procDataFileIDs)


for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    MADataFileID = [procDataFileID(1:end-12) 'MAData.mat'];
    trainingDataFileID = [procDataFileID(1:end-12) 'TrainingData.mat'];
    load(MADataFileID);
    load(trainingDataFileID);
    Diff_MALabels = diff(MALabels);
    Changepoint = find(Diff_MALabels~=0);
%     Changepoint = Changepoint+1;
    for CP = 1:1:length(Changepoint)
        if CP==1
            MASize(CP) = Changepoint(CP);
        else
            MASize(CP) = Changepoint(CP) - Changepoint(CP-1);
        end

            if MALabels(Changepoint(CP)) == 1
                if MASize(CP) < 10
                    MicroAr(CP) = 1; 
                elseif MASize(CP) >=10
                    MicroAr(CP) = 0;
                end           
            elseif MALabels(Changepoint(CP)) == 2
                MicroAr(CP) = 2;
    
            else 
                MicroAr(CP) = 3;
    
            end


    end

end

end
