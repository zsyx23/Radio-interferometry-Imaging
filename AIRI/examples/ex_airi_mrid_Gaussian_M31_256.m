clear 
clc

path = fileparts(mfilename('fullpath'));
cd(path)
cd ..

config = ['.', filesep, 'config', filesep, 'airi_sim_Gaussian_M31_256.json'];
dataFile = ['.', filesep, 'data', filesep, 'data_Gaussian_M31_256_30dB.mat'];
groundtruth = ['.', filesep, 'data', filesep, 'M31_256.fits'];
resultPath = ['.', filesep, 'results_Gaussian_M31_256']; 
algorithm = 'airi';
shelf_pth = ['.', filesep, 'airi_denoisers', filesep, 'shelf_mrid.csv'];
RunID = 1;

run_imager(config, 'dataFile', dataFile, 'algorithm', algorithm, 'resultPath', resultPath, 'dnnShelfPath', shelf_pth, 'groundtruth', groundtruth, 'runID', RunID)