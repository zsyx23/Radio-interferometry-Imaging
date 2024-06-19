clear all
close all
clc

addpath simulated_data/data/
addpath simulated_data/
addpath lib/
addpath lib/irt/
addpath alg

try
    run('./lib/irt/setup.m');
end


generate_groundtruth = 1;
generate_undersampled_cube = 1;

generate_data = 1;
generate_uv = 1;
save_data = 0;
load_data = 0;
save_full_operator = 0;
free_memory = 1;

generate_measurements = 1;

solve_minimization = 1;

compute_Anorm = 1;
compute_Anorm_ch = 0;
generate_eps_nnls = 0;

% choose the solver you want to use
solve_HS = 1; % wide-band solver (rwLRJS) applied on all the channels
solve_1B = 0; % single-channel solver (SARA) applied for each channel separately
load_temp_result = 0;

%% Generate or Load ground-truth
if generate_groundtruth
    file_name = './simulated_data/data/W28_256.fits';
    Nx=256;
    Ny=256;
    c = 60;
    f = linspace(1,2,c); % frequency bandwidth from 1 to 2 GHz with 15 channels
    emission_lines = 1; % insert emission lines on top of the continuous spectra
    [x0,X0] = Generate_cube(file_name,Ny,Nx,f,emission_lines);
elseif load_data
    load('./simulated_data/data/data_60_HI_final3.mat');
end

if generate_undersampled_cube
    % new desired dimensions
    Nx=256;
    Ny=256;
    % undersampling factor for the channels
    unds = 4; % take 1/unds images
    [x0,X0,f,c] = Generate_undersampled_cube(x0,f,Ny,Nx,c,unds);
end

ch = [1 : c]; % number of channels loaded (note that this can be one).
clear c;

%%
percentage = 0.3; [0.3 0.05 0.1 0.02 0.01 0.2]; % sampling rate
for k = 1 : length(percentage)
    
    %diary(['diary_HI_60_' num2str(percentage(k))]);
    
    %% Generate or Load subsampling mask
    if generate_data
        Generate_data
    end
    
    %%
    input_snr = 40; % noise level (on the measurements)
    num_tests = 1;
    % Generate or Load subsampling mask
    if generate_measurements
        Generate_Measurements
        if save_data
            save(['./simulated_data/data/y_60_HI_finalNEW_' num2str(percentage(k)) '.mat'],'-v7.3', 'y0_t', 'y_t', 'y0b_t', 'yb_t', 'epsilonT', 'epsilonT1', 'sigma_noise');
        end
    elseif load_data
        load(['y_60_HI_finalNEW_' num2str(percentage(k)) '.mat']);
        disp('data loaded successfully ;)')
    end
    
    %%
    if solve_minimization
        
        Solver_simulated_data
    end
    
end
