clear;
close all;
clc

path = fileparts(mfilename('fullpath'));
cd(path)
cd ..

average_times= 1 ;

per=0.1:0.1:0.9;
SNR_sum_AIRI=[]; 
SNRlog_sum_AIRI=[];  
TIME_sum_AIRI=[];

SNR_sum_ave_AIRI=zeros(1,length(per)); 
SNRlog_sum_ave_AIRI=zeros(1,length(per)); 
TIME_sum_ave_AIRI=zeros(1,length(per));

      dataFile = ['.', filesep, 'data', filesep, 'AIRI_Gaussian_M31_256_25dB_40.mat'];   
      config = ['.', filesep, 'config', filesep, 'airi_sim_Gaussian_M31_256_per_40.json'];
   

groundtruth = ['.', filesep, 'data', filesep, 'M31_256.fits'];
resultPath = ['.', filesep, 'results_Gaussian_M31_256_per']; 
algorithm = 'airi';
shelf_pth = ['.', filesep, 'airi_denoisers', filesep, 'shelf_mrid.csv'];
RunID = 1;

% run_imager_per(config, 'dataFile', dataFile, 'algorithm', algorithm, 'resultPath', resultPath, 'dnnShelfPath', shelf_pth, 'groundtruth', groundtruth, 'runID', RunID)
tic
[MODEL]=run_imager_per(config, 'dataFile', dataFile, 'algorithm', algorithm, 'resultPath', resultPath, 'dnnShelfPath', shelf_pth, 'groundtruth', groundtruth, 'runID', RunID);
TIME_AIRI=toc;



