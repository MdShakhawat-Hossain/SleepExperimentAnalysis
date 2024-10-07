
clear; clc; close all;
%% Animal IDs
FP_animalIDs =    {'NEACh007','NEACh008'};
%% make sure the code repository and data are present in the current directory
firstHrs = "false";
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
%%
if exist('GLMDataMatrix.mat','file') == 2
    load('GLMDataMatrix.mat')
else
    GLMDataMatrix = [];
end
%%
runFromStart = 'y';
for kk = 1:length(FP_animalIDs)
    if isfield(GLMDataMatrix,(FP_animalIDs{1,kk})) == false || isfield(GLMDataMatrix.(FP_animalIDs{1,kk}),'GLMModelData') == false || strcmp(runFromStart,'y') == true
          [GLMDataMatrix] = AnalyzeGLM_Data(FP_animalIDs{1,kk},rootFolder,GLMDataMatrix);
    end
    multiWaitbar('Separating Data for GLM','Value',kk/length(FP_animalIDs));
end