function [] = MainScript_FP_GRABNE_firstHrs()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
% Purpose: Generates KLT's main and supplemental figs for Turner et al. eLife2020
%
% Scripts used to pre-process the original data are located in the folder "Pre-Processing Scripts".
% Functions that are used in both the analysis and pre-processing are located in the analysis folder.
%________________________________________________________________________________________________________________________
% addpath(genpath('C:\Users\mfh5734\OneDrive - The Pennsylvania State University\Documents\Research_Codes\FiberPhotometry\Data-Analysis-master'))
clear; clc; close all;

%% Animal ID
FP_animalIDs =  {'NEACh007','NEACh008'}; %
%% make sure the code repository and data are present in the current directory
% currentFolder = 'H:\Sleep_GCaMP7s_ChATCre\';
% these data are from first hours
firstHrs = 'true';
%
currentFolder = pwd;
addpath(genpath(currentFolder));
fileparts = strsplit(currentFolder,filesep);
if ismac
    rootFolder = fullfile(filesep,fileparts{1:end});
    delim = '/';
else
    rootFolder = fullfile(fileparts{1:end});
    delim = '\';
end
% add root folder to Matlab's working directory
addpath(genpath(rootFolder))
%% run the data analysis. The progress bars will show the analysis progress
rerunAnalysis = 'y';
saveFigs = 'y';
if exist('AnalysisResults_firstHrs.mat','file') ~= 2 || strcmp(rerunAnalysis,'y') == true
%     multiWaitbar('Analyzing sleep probability',0,'Color','B'); pause(0.25);
    % run analysis and output a structure containing all the analyzed data
    [AnalysisResults_firstHrs] = AnalyzeData(rootFolder,FP_animalIDs);
    multiWaitbar('CloseAll');
else
    disp('Loading analysis results and generating figures...'); disp(' ')
    load('AnalysisResults_firstHrs.mat')
end
%% supplemental figure panels
% [AnalysisResults_firstHrs] = Fig8_S1(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig7_S3(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig7_S2(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig7_S1(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
%% [AnalysisResults_firstHrs] = Fig6_S1(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig5_S1(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig3_S5(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig3_S4(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig3_S3(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig3_S2(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig3_S1(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig2_S2(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig2_S1(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig1_S9(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig1_S8(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig1_S7(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig1_S6(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig1_S5(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig4_S1_GRABNE_Rhodamine(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig4_S1_GCaMP7s(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);

% [AnalysisResults_firstHrs] = Fig4_FP_Bilateral_No_REM_TRITCChAT_movement(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig4_FP_Bilateral_No_REM_TRITCChAT_movement_compare(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);


% [AnalysisResults_firstHrs] = Fig1_S1(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
%% main figure panels
% [AnalysisResults_firstHrs] = Fig1_S4_FP_Stats_GRABNE(rootFolder,saveFigs,delim,AnalysisResults_firstHrs,firstHrs);
% [AnalysisResults_firstHrs] = Fig1_S3_Whisk_GRABNE(rootFolder,saveFigs,delim,AnalysisResults_firstHrs,firstHrs);
[AnalysisResults_firstHrs] = Fig1_S2_Stim_GRABNE_Response(rootFolder,saveFigs,delim,AnalysisResults_firstHrs,firstHrs,FP_animalIDs);
% [AnalysisResults_firstHrs] = Fig1_S2_Stim_GRABNE(rootFolder,saveFigs,delim,AnalysisResults_firstHrs,firstHrs,FP_animalIDs);
% [AnalysisResults_firstHrs] = Fig1_S2_Stim_GRABNE_SingleAnimal(rootFolder,saveFigs,delim,AnalysisResults_firstHrs,firstHrs,FP_animalIDs);
% [AnalysisResults_firstHrs] = Fig1_S2_Stim_GRABNE_compare(rootFolder,saveFigs,delim,AnalysisResults_firstHrs,firstHrs);
% [AnalysisResults_firstHrs] = PlotCOherence_GRABNE(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = PlotpowerSpectrum_GRABNE(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = PlotCrossCorrelation_GRABNE(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig4_FP_Transition_GRABNE(rootFolder,saveFigs,delim,AnalysisResults_firstHrs,firstHrs);
% [AnalysisResults_firstHrs] = Fig1_S5_FP_Stats_GRABNE(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig1_S6_FP_Stats_GRABNE(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
% [AnalysisResults_firstHrs] = Fig1_S7_FP_Stats_GRABNE(rootFolder,saveFigs,delim,AnalysisResults_firstHrs);
%% fin.
disp('MainScript Analysis - Complete'); disp(' ')
end

function [AnalysisResults_firstHrs] = AnalyzeData(rootFolder,FP_animalIDs)
% FP animal IDs
% FP_animalIDs = {'NEACh003'};

saveFigs = 'y';
if exist('AnalysisResults_firstHrs.mat','file') == 2
    load('AnalysisResults_firstHrs.mat')
else
    AnalysisResults_firstHrs = [];
end
% these data are from first hours
firstHrs = 'true';
%% Block [1] Analyze the arousal-state probability of trial duration and resting events (IOS)
% runFromStart = 'n';
% for aa = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,aa})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,aa}),'SleepProbability') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults_firstHrs] = AnalyzeAwakeProbability_GRBANE(FP_animalIDs{1,aa},saveFigs,rootFolder,AnalysisResults_firstHrs,firstHrs);
%     end
%     multiWaitbar('Analyzing sleep probability','Value',aa/length(FP_animalIDs));
% end
%% Block [2] Analyze the arousal-state distribution of different behavioral measurements (IOS)
% runFromStart = 'n';
% for bb = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,bb})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,bb}),'BehaviorDistributions') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults_firstHrs] = AnalyzeBehavioralDistributions_GRABNE(FP_animalIDs{1,bb},rootFolder,AnalysisResults_firstHrs,firstHrs);
%     end
%     multiWaitbar('Analyzing behavioral distributions','Value',bb/length(FP_animalIDs));
% end
%% Block [3] Analyze the transitions between different arousal-states (IOS)
% runFromStart = 'n';
% for dd = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,dd})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,dd}),'Transitions') == false || strcmp(runFromStart,'y') == true
%           [AnalysisResults_firstHrs] = AnalyzeTransitionalAverages_FP_Movements_GRABNE(FP_animalIDs{1,dd},saveFigs,rootFolder,AnalysisResults_firstHrs,firstHrs);
%     end
%     multiWaitbar('Analyzing behavioral transitions triggered changes','Value',dd/length(FP_animalIDs));
% end
%% Block [4] Analyze the hemodynamic signal [HbT] during different arousal states (IOS)
% runFromStart = 'y';
% for ff = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,ff})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,ff}),'MeanRhodamine') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults_firstHrs] = AnalyzeMeanRhodamine_FP_GRABNE(FP_animalIDs{1,ff},rootFolder,AnalysisResults_firstHrs,firstHrs);
%     end
%     multiWaitbar('Analyzing behavioral Rhodamine','Value',ff/length(FP_animalIDs));
% end
%% Block [5] Analyze the hemodynamic signal [GCaMP7s] during different arousal states (IOS)
% runFromStart = 'y';
% for ff = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,ff})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,ff}),'MeanGFP') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults_firstHrs] = AnalyzeMeanGFP_FP_GRABNE(FP_animalIDs{1,ff},rootFolder,AnalysisResults_firstHrs,firstHrs);
%     end
%     multiWaitbar('Analyzing behavioral GFP','Value',ff/length(FP_animalIDs));
% end
%% Block [6] Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
runFromStart = 'y';
for pp = 1:length(FP_animalIDs)
    if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,pp})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,pp}),'Whisk') == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,pp}),'Stim') == false ||strcmp(runFromStart,'y') == true
        [AnalysisResults_firstHrs] = AnalyzeEvokedResponses_GRABNE(FP_animalIDs{1,pp},rootFolder,AnalysisResults_firstHrs,firstHrs);
    end
    multiWaitbar('Analyzing evoked responses','Value',pp/length(FP_animalIDs));
end
%% Block [7] Analyze the spectral coherence between bilateral hemodynamic [HbT] and neural signals (IOS)
% runFromStart = 'n';
% for jj = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,jj})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,jj}),'Coherence') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults_firstHrs] = AnalyzeCoherence_Pupil_GRABNE(FP_animalIDs{1,jj},rootFolder,AnalysisResults_firstHrs);
%     end
%     multiWaitbar('Analyzing coherence','Value',jj/length(FP_animalIDs));
% end
%% Block [8] Analyze the spectral coherence between neural-hemodynamic [HbT] signals (IOS)
% runFromStart = 'n';
% for jj = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,jj})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,jj}),'NeuralHemoCoherence') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults_firstHrs] = AnalyzeNeuralHemoCoherence(FP_animalIDs{1,jj},saveFigs,rootFolder,AnalysisResults_firstHrs);
%     end
%     multiWaitbar('Analyzing neural-hemo coherence','Value',jj/length(FP_animalIDs));
% end
%% Block [9] Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
% runFromStart = 'n';
% for kk = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,kk})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,kk}),'PowerSpectrum') == false || strcmp(runFromStart,'y') == true
% %         [AnalysisResults_firstHrs] = AnalyzePowerSpectrum(FP_animalIDs{1,kk},rootFolder,AnalysisResults_firstHrs);
%           [AnalysisResults_firstHrs] = AnalyzePowerSpectrum_Pupil_GRABNE(FP_animalIDs{1,kk},rootFolder,AnalysisResults_firstHrs);
%     end
%     multiWaitbar('Analyzing power spectra','Value',kk/length(FP_animalIDs));
% end
%% Block [10] Analyze Pearson's correlation coefficient between bilateral hemodynamic [HbT] and neural signals (IOS)
% runFromStart = 'y';
% for mm = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,mm})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,mm}),'CorrCoeff') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults_firstHrs] = AnalyzeCorrCoeffs(FP_animalIDs{1,mm},rootFolder,AnalysisResults_firstHrs);
%     end
%     multiWaitbar('Analyzing Pearson''s correlation coefficients','Value',mm/length(FP_animalIDs));
% end
%% Block [11] Analyze the cross-correlation between neural activity and hemodynamics [HbT] 
% runFromStart = 'n';
% for nn = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,nn})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,nn}),'CrossCorrelation') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults_firstHrs] = AnalyzeCrossCorrelation_GRABNE(FP_animalIDs{1,nn},rootFolder,AnalysisResults_firstHrs);
%     end
%     multiWaitbar('Analyzing cross correlation','Value',nn/length(FP_animalIDs));
% end

%% Block [14] Analyze the relationship between gamma-band power and hemodynamics [HbT] (IOS)
% runFromStart = 'y';
% for qq = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,qq})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,qq}),'TRITCvsGamma') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults_firstHrs] = AnalyzeCBVGammaRelationship(FP_animalIDs{1,qq},rootFolder,AnalysisResults_firstHrs);
%     end
%     multiWaitbar('Analyzing CBV-Gamma relationship','Value',qq/length(FP_animalIDs));
% end
 %% Block [15] Analyze the relationship between gamma-band power and GCaMP7s
% runFromStart = 'y';
% for qq = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,qq})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,qq}),'GCaMP7svsGamma') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults_firstHrs] = AnalyzeGCaMP7sGammaRelationship(FP_animalIDs{1,qq},rootFolder,AnalysisResults_firstHrs);
%     end
%     multiWaitbar('Analyzing GCaMP7s-Gamma relationship','Value',qq/length(FP_animalIDs));
% end
%% Block [16] Analyze the relationship between HBT and GCaMP7s
% runFromStart = 'y';
% for qq = 1:length(FP_animalIDs)
%     if isfield(AnalysisResults_firstHrs,(FP_animalIDs{1,qq})) == false || isfield(AnalysisResults_firstHrs.(FP_animalIDs{1,qq}),'TRITCvsGCaMP7s') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults_firstHrs] = AnalyzeHbTGCaMP7sRelationship(FP_animalIDs{1,qq},rootFolder,AnalysisResults_firstHrs);
%     end
%     multiWaitbar('Analyzing CBV-GCaMP7s relationship','Value',qq/length(FP_animalIDs));
% end
%% fin.
disp('Loading analysis results and generating figures...'); disp(' ')

end
