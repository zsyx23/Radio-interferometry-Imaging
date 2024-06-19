clear all
close all
clc

addpath real_data/data/
addpath real_data/
addpath lib/
addpath lib/irt/
addpath alg

try
    run('./lib/irt/setup.m');
end

compute_Anorm = 1;
compute_Anorm_ch = 0;
generate_measurement_operator = 1;
generate_eps_nnls = 1;

save_data = 0;
load_data = 0;
save_full_operator = 0;
free_memory = 1;

solve_minimization = 1;

% choose the solver you want to use
solve_HS = 1; % wide-band solver (rwLRJS) applied on all the channels
solve_1B = 0; % single-channel solver (SARA) applied for each channel separately
load_temp_result = 0;

%% Config parameters if uploading saved measurement operator data
if compute_Anorm || ~generate_measurement_operator
    Nx = 1024;
    Ny = 512;
    N = Nx * Ny;
    ox = 2; % oversampling factors for nufft
    oy = 2; % oversampling factors for nufft
    Kx = 8; % number of neighbours for nufft
    Ky = 8; % number of neighbours for nufft
end

%% Generate or Load real data
if generate_measurement_operator
    measurement_operator_real_data;
else
    if load_data
        load('./real_data/data/CYG_data.mat');
    end
    ch = [1 : size(yb,2)]
    [A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);
end

%%
if solve_minimization
    
    Solver_real_data;
    
end

