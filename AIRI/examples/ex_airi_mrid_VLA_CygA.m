clear 
clc
close all;

path = fileparts(mfilename('fullpath'));
cd(path)
cd ..

config = ['.', filesep, 'config', filesep, 'airi_sim_VLA.json'];
dataFile = ['.', filesep, 'data', filesep, 'data_VLA_25dB_CygA.mat'];
groundtruth = ['.', filesep, 'data', filesep, 'data_CygA_1024.fits'];
resultPath = ['.', filesep, 'results_VLA']; 
algorithm = 'airi'; 
shelf_pth = ['.', filesep, 'airi_denoisers', filesep, 'shelf_mrid.csv'];
RunID = 1;

tic
[MODEL]=run_imager_per(config, 'dataFile', dataFile, 'algorithm', algorithm, 'resultPath', resultPath, 'dnnShelfPath', shelf_pth, 'groundtruth', groundtruth, 'runID', RunID);
TIME_AIRI=toc




