function PatchPupilArea_FP(rawDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: The pupil camera occasionally drops packets of frames. We can calculate the difference in the number
%          of expected frames as well as the indeces that LabVIEW found the packets lost.
%________________________________________________________________________________________________________________________
    for qq = 1:size(rawDataFileIDs,1)
        rawDataFileID = rawDataFileIDs(qq,:);
        load(rawDataFileID)
        disp(['Patching missing pupil data of file ' num2str(qq) '/' num2str(size(rawDataFileIDs,1))]); disp(' ')
        if isfield(RawData.data.Pupil,'pupilPatch') == false
            [animalID,~,fileID] = GetFileInfo_JNeurosci2022(rawDataFileIDs);
            % expected number of frames based on trial duration and sampling rate
            expectedSamples = RawData.notes.trialDuration_long*RawData.notes.pupilCamSamplingRate;
            droppedFrameIndex = RawData.notes.droppedPupilCamFrameIndex;
            sampleDiff = expectedSamples - length(RawData.data.Pupil.pupilArea);
            framesPerIndex = ceil(sampleDiff/length(droppedFrameIndex));
            blinks = RawData.data.Pupil.blinkInds;
            %% patch NaN values
            pupilArea = RawData.data.Pupil.pupilArea;
            nanLogical = isnan(pupilArea);
            nanIndex = find(nanLogical == 1);
            if sum(nanLogical) > 1 && sum(nanLogical) < 1000
                while sum(nanLogical) >= 1
                    pupilArea = fillmissing(pupilArea,'movmedian',3);
                    nanLogical = isnan(pupilArea);
                end
                % check length of missing data. If there's periods > than 1 second of continuous NaN then mark the file as bad
                testPatch = fillmissing(RawData.data.Pupil.pupilArea,'movmedian',31);
                if sum(isnan(testPatch)) > 1
                    RawData.data.Pupil.diameterCheck = 'n';
                    RawData.data.Pupil.diameterCheckComplete = 'y';
                else
                    nanFigure = figure;
                    plot((1:length(pupilArea))./RawData.notes.pupilCamSamplingRate,pupilArea,'k');
                    hold on
                    for aa = 1:length(nanIndex)
                        x1 = xline(nanIndex(1,aa)/RawData.notes.pupilCamSamplingRate,'r');
                    end
                    if isempty(droppedFrameIndex) == false
                        for bb = 1:length(droppedFrameIndex)
                            x2 = xline(droppedFrameIndex(1,bb)/RawData.notes.pupilCamSamplingRate,'g');
                        end
                        legend([x1,x2],'NaNs','Dropped frames')
                    end
                    title([animalID ' ' strrep(fileID,'_',' ')])
                    legend(x1,'NaNs')
                    xlabel('Time (sec)')
                    ylabel('Pupil area (pixels')
                    axis tight
                    % save the file to directory.
                    [pathstr,~,~] = fileparts(cd);
                    dirpath = [pathstr '/Figures/Pupil Data Patching/'];
                    if ~exist(dirpath,'dir')
                        mkdir(dirpath);
                    end
                    savefig(nanFigure,[dirpath animalID '_' fileID '_PupilPatch'])
                    close(nanFigure)
                end
            end
            %% patch missing frames now that NaN are gone
            if ~isempty(droppedFrameIndex)
                addedFrames = 0;
                % each dropped index
                for cc = 1:length(droppedFrameIndex)
                    % for the first event, it's okay to start at the actual index
                    if cc == 1
                        leftEdge = droppedFrameIndex(1,cc);
                    else
                        % for all other dropped frames after the first, we need to correct for the fact that index is shifted right.
                        leftEdge = droppedFrameIndex(1,cc) + ((cc - 1)*framesPerIndex);
                    end
                    % set the edges for the interpolation points. we want n number of samples between the two points,vthe left and
                    % right edge values. This equates to having a 1/(dropped frames + 1) step size between the edges.
                    rightEdge = leftEdge + 1;
                    patchFrameInds = leftEdge:(1/(framesPerIndex + 1)):rightEdge;
                    % concatenate the original data for the first index, then the new patched data for all subsequent
                    % indeces. Take the values from 1:left edge, add in the new frames, then right edge to end.
                    if cc == 1
                        patchFrameVals = interp1(1:length(pupilArea),pupilArea,patchFrameInds); % linear interp
                        snipPatchFrameVals = patchFrameVals(2:end - 1);
                        try
                            patchedPupilArea = horzcat(pupilArea(1:leftEdge),snipPatchFrameVals,pupilArea(rightEdge:end));
                        catch
                            patchedPupilArea = horzcat(pupilArea(1:end),snipPatchFrameVals);
                        end
                    else
                        patchFrameVals = interp1(1:length(patchedPupilArea),patchedPupilArea,patchFrameInds); % linear interp
                        snipPatchFrameVals = patchFrameVals(2:end - 1);
                        patchedPupilArea = horzcat(patchedPupilArea(1:leftEdge),snipPatchFrameVals,patchedPupilArea(rightEdge:end));
                    end
                    addedFrames = addedFrames + length(snipPatchFrameVals);
                    addedFrameIndex(1,cc) = addedFrames;
                    if cc == 1
                        for qq = 1:length(blinks)
                            if leftEdge < blinks(1,qq)
                                shiftedBlinks(1,qq) = blinks(1,qq) + addedFrames;
                            else
                                shiftedBlinks(1,qq) = blinks(1,qq);
                            end
                        end
                    else
                        for qq = 1:length(blinks)
                            if leftEdge < blinks(1,qq)
                                shiftedBlinks(1,qq) = shiftedBlinks(1,qq) + addedFrames;
                            else
                                shiftedBlinks(1,qq) = shiftedBlinks(1,qq);
                            end
                        end
                    end
                end
                % due to rounding up on the number of dropped frames per index, we have a few extra frames. Snip them off.
                patchedPupilArea = patchedPupilArea(1:expectedSamples);
                if isempty(blinks) == false
                    RawData.data.Pupil.shiftedBlinks = shiftedBlinks;
                    RawData.data.Pupil.addedFrameIndex = addedFrameIndex;
                end
            else
                patchedPupilArea = pupilArea(1:expectedSamples);
            end
            RawData.data.Pupil.pupilArea = patchedPupilArea;
            RawData.data.Pupil.pupilPatch = 'y';
            save(rawDataFileID,'RawData','-v7.3')
        end
    end

end
