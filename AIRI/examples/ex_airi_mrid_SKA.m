clear 
clc
close all;

path = fileparts(mfilename('fullpath'));
cd(path)
cd ..

config = ['.', filesep, 'config', filesep, 'airi_sim_SKA.json'];
dataFile = ['.', filesep, 'data', filesep, 'data_SKA_25dB.mat'];
groundtruth = ['.', filesep, 'data', filesep, 'W28.fits'];
resultPath = ['.', filesep, 'results_SKA']; 
algorithm = 'airi'; 
shelf_pth = ['.', filesep, 'airi_denoisers', filesep, 'shelf_mrid.csv'];
RunID = 1;

tic
[MODEL]=run_imager_per(config, 'dataFile', dataFile, 'algorithm', algorithm, 'resultPath', resultPath, 'dnnShelfPath', shelf_pth, 'groundtruth', groundtruth, 'runID', RunID);
TIME_AIRI=toc



