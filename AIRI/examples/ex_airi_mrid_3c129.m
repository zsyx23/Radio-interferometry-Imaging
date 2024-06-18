clear 
clc

path = fileparts(mfilename('fullpath'));
cd(path)
cd ..

config = ['.', filesep, 'config', filesep, 'airi_sim_3c129.json'];
dataFile = ['.', filesep, 'data', filesep, 'data_3c129.mat'];
groundtruth = ['.', filesep, 'data', filesep, 'W28.fits'];
resultPath = ['.', filesep, 'results_3c129']; 
algorithm = 'airi'; 
shelf_pth = ['.', filesep, 'airi_denoisers', filesep, 'shelf_mrid.csv'];
RunID = 1;

run_imager(config, 'dataFile', dataFile, 'algorithm', algorithm, 'resultPath', resultPath, 'dnnShelfPath', shelf_pth, 'groundtruth', groundtruth, 'runID', RunID)


