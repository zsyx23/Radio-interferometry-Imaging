clear;
close all;

load('data_test_main_all.mat');

addpath data/
addpath data/img
addpath data/uvw
addpath data/vis
addpath lib/
addpath alg/

%%%%  The program runs up to here,
%%%%  and then manually adds 
%%%%  CubeHelix-master, fessler, function and iuwt paths.

fprintf('Generating new data ... \n\n');
%% image loading
% [im, N, Ny, Nx] = util_read_image(image_file_name);

importdata M31_256.fits ;
im=double(ans);
im=im/(max(max(im)));
[Nx,Ny]=size(im);
im=imresize(im,256./Nx);
[Nx,Ny]=size(im);N=Nx*Ny;

figure; imagesc(real(log10((im))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-2 0]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

%% 
input_snr = 25 ; % noise level (on the measurements)
param_sampling.p=0.4;
Maxiter= 200 ;soft_threshold_values=1e-6; forward_descent_step_size=0.1 ;
soft_threshold_values_admm=1e-6;
lambda_IUWT= 1e-5 ;  gamma_IUWT = 0.5 ;  iter_max_IUWT  =  200 ;
gamma_WTF= 0.1 ;  iter_max_WTF= 100 ;
f_log_min=-2;
f_log_max=0;
f_res=0.005;
f_err=0.05;

%% generate the sampling pattern
param_sampling.N = N; % number of pixels in the image
param_sampling.Nox = ox*Nx; % number of pixels in the image
param_sampling.Noy = oy*Ny; % number of pixels in the image

[uw, vw, ~] = util_gen_sampling_pattern(sampling_pattern, param_sampling);


if use_symmetric_fourier_sampling
    uw = [uw; -uw];
    vw = [vw; -vw];
end

%% 
figure;plot(uw,vw,'.k');
axis equal;
axis([-pi pi -pi pi]);
xlabel('normalised frequency','FontName','Times New Roman','FontSize',16,'LineWidth',2);% 字体大小30
ylabel('normalised frequency','FontName','Times New Roman','FontSize',16,'LineWidth',2);
set(gca,'position',[0.08,0.1,0.88,0.88])
set(gcf,'unit','normalized','position',[0.2,0.2,0.3,0.45]);

%% compute weights
param_precond.N = N; % number of pixels in the image
param_precond.Nox = ox*Nx; % number of pixels in the image
param_precond.Noy = oy*Ny; % number of pixels in the image
[aWw] = util_gen_preconditioning_matrix(uw, vw, param_precond);

%% set the blocks structure
nWw = ones(length(uw), 1);
[u, v, ~, uvidx, aW, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);


%% measurement operator initialization

fprintf('Initializing the NUFFT operator\n\n');
tstart = tic;
[A, At, G, W, Gw] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
tend = toc(tstart);

fprintf('Initialization runtime: %ds\n\n', ceil(tend));
R = length(v);

y0 = cell(num_tests, 1);
y = cell(num_tests, 1);
aY = cell(num_tests, 1);
epsilon = cell(num_tests, 1);
epsilons = cell(num_tests, 1);
epsilonT = cell(num_tests, 1);
epsilonTs = cell(num_tests, 1);
y0f = cell(num_tests, 1);
yf = cell(num_tests, 1);
input_snr_v = cell(num_tests, 1);

%% generate noisy input data
for k = 1:num_tests
    [y0{k}, y{k}, y0f{k}, yf{k}, aY{k}, input_snr_v{k}] = util_gen_input_data(im, G, W, A, input_snr, use_different_per_block_input_snr, per_block_input_snr_delta, uvidx);
    
    if use_symmetric_fourier_sampling
        y0f{k} = [y0f{k}(uvidx{k}(1:end/2)); conj(y0f{k}(uvidx{k}(1:end/2)))];
        yf{k} = [yf{k}(uvidx{k}(1:end/2)); conj(yf{k}(uvidx{k}(1:end/2)))];
        for j = 1:R
            y{k}{j} = [y{k}{j}(uvidx{k}(1:end/2)); conj(y{k}{j}(uvidx{k}(1:end/2)))];
            y0{k}{j} = [y0{k}{j}(uvidx{k}(1:end/2)); conj(y0{k}{j}(uvidx{k}(1:end/2)))];
            aY{k}{j} = [aY{k}{j}(uvidx{k}(1:end/2)); conj(aY{k}{j}(uvidx{k}(1:end/2)))];
        end
    end
    
    
    [epsilonT{k}, epsilonTs{k}, epsilon{k}, epsilons{k}] = util_gen_L2_bounds(y{k}, ...
        input_snr, [], l2_ball_definition, stopping_criterion, use_same_stop_criterion, param_l2_ball);
end

%% gridding
if use_gridded_data
    [yT, T, Tw] = util_compute_gridding_data(y, G, Gw);
else
    T = G;
    Tw = Gw;
    yT = y;
end

%% sparsity operator definition

[Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
[Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);


%% compute the operator norm
fprintf('Computing operator norms ...\n');

fprintf('Natural W ...\n');
evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-6, 200, verbosity);

fprintf('No natural W ...\n');
Gw_ = spdiags(1./cell2mat(nW), 0, length(nWw), length(nWw)) * Gw;
Gwt_ = Gw_';
evl_no_natw = op_norm(@(x) Gw_ * A(x), @(x) At(Gwt_ * x), [Ny, Nx], 1e-6, 200, verbosity);
clear Gw_ Gwt_;

fprintf('Preconditioning ...\n');
evl_precond = op_norm(@(x) sqrt(cell2mat(aW)) .* (Gw * A(x)), @(x) At(Gw' * (sqrt(cell2mat(aW)) .* x)), [Ny, Nx], 1e-6, 200, verbosity);

evl_blocks = zeros(R, 1);

if compute_block_op_norm
    % maximum eigenvalue of operator
    for q = 1:R
        No = size(W{1}, 1);
        Tw_ = spalloc(size(T{q}, 1), No, size(T{q}, 2) * 16);
        Tw_(:, W{q}) = T{q};
        fprintf('\nComputing operator norm: block %i \n', q)
        evl_blocks(q) = op_norm(@(x) sqrt(aW{q}) .* (Tw_ * A(x)), @(x) At(Tw_' * (sqrt(aW{q}) .* x)), [Ny, Nx], 1e-6, 200, verbosity);
        clear Tw_;
    end
end

%% save data
if save_data_on_disk == 1
    fprintf('Saving new data ... \n');
    
    if save_data_on_disk
        file_start = save_dataset_number;
        while exist(sprintf('%s%s_input_data.mat', save_path, int2str(file_start)), 'file') || ...
                exist(sprintf('$s%s_input_data_config.mat', save_path, int2str(file_start)), 'file')
            file_start = file_start + 1;
        end
        if file_start ~= save_dataset_number;
            fprintf('WARNING: Saving new data in file %d instead of %d \n\n', file_start, save_dataset_number);
        end
        
%         save(sprintf('%s%s_input_data', save_path, int2str(file_start)), '-v7.3', ... % 'G', 'Gw', 'W'
%             'y0f', 'yf', 'y', 'y0', 'nWw'); % do not save G, it is large and is faster to compute than saving
%         save(sprintf('%s%s_input_data_config', save_path, int2str(file_start)), 'N', 'Ny', ...
%             'Nx', 'uw', 'vw', 'u', 'v', 'uvidx', 'im', 'sampling_pattern', 'param_sampling', 'input_snr', 'input_snr_v', 'image_file_name', 'num_tests', 'use_real_visibilities');
        
        if strcmp(sampling_pattern, 'file')
            for k = 1:num_tests
                vis = yf{k};
%                 save(sprintf('%s%s_vis_data_t%s', save_path, int2str(file_start), int2str(k)), 'vis');
                clear vis;
                
                vis = yf{k} - y0f{k};
%                 save(sprintf('%s%s_vis_data_noise_t%s', save_path, int2str(file_start), int2str(k)), 'vis');
                clear vis;
            end
        end
    end
end


if gen_data == 2
    fprintf('Using data from workspace ... \n\n');
end

if free_memory
    % free memory of the whose Gw or the split G are not needed
    try
        if ~run_admm_bpcon && ~run_sdmm_bpcon && ~run_fb_nnls && ~run_krylov_nnls
            clear Gw;
            clear Tw;
        end
        
        if ~run_pdfb_bpcon_par_sim_rescaled_precond_wave_par && ~run_pdfb_bpcon_par_sim && ~run_pdfb_bpcon_par_sim_rescaled && ~run_pdfb_bpcon_par_sim_rescaled_precond &&...
                ~run_pdfb_bpcon_par_sim_rand_rescaled && ~run_pdfb_bpcon_par_sim_block_rand_rescaled && ...
                ~run_pdfb_bpcon_dist && ~run_pdfb_bpcon_dist_rescaled && ~run_admm_bpconpar && ~run_sdmm_bpconpar && ...
                ~run_pdfb_bpcon_par_sim_rescaled_rec_async && ...
                ~run_pdfb_bpcon_par_sim_rand_rescaled_nonuniform_p && ...
                ~run_admm_bpconpar_wavepar
            clear G;
            clear T;
            clear W;
        end
    end
end



%% compute the solution

% snr = 0;
% sp = 0;
% soln = zeros(Ny, Nx);
% residualn = zeros(Ny, Nx);
% dirtyn = zeros(Ny, Nx);
% asol = result_st.sol{i};
% ay = y{i};
% aL1_v = result_st.L1_v{i};
% aL2_v = result_st.L2_v{i};
% aL1_vp = result_st.L1_vp{i};
% aL2_vp = result_st.L2_vp{i};
% % ay0 = y0{i};
% aepsilon = epsilon{i};
% aepsilonT = epsilonT{i};
% R = length(ay);
% ys = A(asol);
% 
% asol = result_st.sol{i};
% residual = At(Gw' * (cell2mat(ay) - Gw * ys));
% residualn_ = real(residual)/evl;
% residualn = residualn + residualn_;
% 
% dirty = At(Gw' * cell2mat(ay));
% dirty = 2 * real(dirty);

%% PD
% run_pdfb_bpcon_par_sim_rescaled

param_pdfb.im = im; % original image, used to compute the SNR
param_pdfb.verbose = verbosity; % print log or not
param_pdfb.nu1 = 1; % bound on the norm of the operator Psi
param_pdfb.nu2 = evl; % bound on the norm of the operator A*G
param_pdfb.gamma = 1e-6; % convergence parameter L1 (soft th parameter)
param_pdfb.tau = 0.49; % forward descent step size
param_pdfb.rel_obj = 0;1e-3; % stopping criterion
% param_pdfb.max_iter = 20; % max number of iterations
param_pdfb.max_iter = 1 ;
param_pdfb.lambda0 = 1; % relaxation step for primal update
param_pdfb.lambda1 = 1; % relaxation step for L1 dual update
param_pdfb.lambda2 = 1; % relaxation step for L2 dual update
param_pdfb.sol_steps = [inf]; % saves images at the given iterations

param_pdfb.use_reweight_steps = 0;
param_pdfb.use_reweight_eps = 0;
param_pdfb.reweight_steps = [200 inf];
param_pdfb.reweight_rel_obj = 1e-5; % criterion for performing reweighting
param_pdfb.reweight_min_steps_rel_obj = 100;
param_pdfb.reweight_alpha = 0.01;
param_pdfb.reweight_alpha_ff = 0.5;

result_st = [];
result_st.sol = cell(num_tests, 1);
result_st.L1_v = cell(num_tests, 1);
result_st.L1_vp = cell(num_tests, 1);
result_st.L2_v = cell(num_tests, 1);
result_st.L2_vp = cell(num_tests, 1);
result_st.time = cell(num_tests, 1);
result_st.delta_v = cell(num_tests, 1);
result_st.sol_v = cell(num_tests, 1);
result_st.sol_reweight_v = cell(num_tests, 1);
result_st.snr_v = cell(num_tests, 1);

result_st.snr = cell(num_tests, 1);
result_st.sparsity = cell(num_tests, 1);
result_st.no_itr = cell(num_tests, 1);


for i = 1:num_tests
    % wavelet mode is a global variable which does not get transfered
    % to the workes; we need to set it manually for each worker
    dwtmode('per');
    
    fprintf('Test run %i:\n', i);
    
    tstart_a = tic;
    fprintf(' Running pdfb_bpcon_par_sim_rescaled\n');
    [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
        result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, ~, ~, result_st.sol_reweight_v{i}] ...
        = pdfb_bpcon_par_sim_rescaled(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, W, Psi, Psit, Psiw, Psitw, param_pdfb);
    tend = toc(tstart_a);
    fprintf(' pdfb_bpcon_par_sim_rescaled runtime: %ds\n\n', ceil(tend));
    
    result_st.time{i} = tend;
    if ~use_real_visibilities
        error = im - result_st.sol{i};
        result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
    end
    result_st.no_itr{i} = length(result_st.L1_v{i});
    
    wcoef = [];
    for q = 1:length(Psit)
        wcoef = [wcoef; Psit{q}(result_st.sol{i})];
    end
    result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
end

% results_prefix = 'pdfb_bpcon_par_sim_rescaled';
% param_structure_name = 'param_pdfb';

snr = 0;
sp = 0;
soln = zeros(Ny, Nx);
residualn = zeros(Ny, Nx);
dirtyn = zeros(Ny, Nx);
asol = result_st.sol{i};
ay = y{i};
aL1_v = result_st.L1_v{i};
aL2_v = result_st.L2_v{i};
aL1_vp = result_st.L1_vp{i};
aL2_vp = result_st.L2_vp{i};
% ay0 = y0{i};
aepsilon = epsilon{i};
aepsilonT = epsilonT{i};
R = length(ay);
ys = A(asol);

% results
%     script_gen_figures;

% save data
%     script_save_result_data;


asol = result_st.sol{i};
residual = At(Gw' * (cell2mat(ay) - Gw * ys));
residualn_PD = real(residual)/evl;
% residualn = residualn + residualn_;

% figure;imagesc(log10(max(im, 1e-8))); colorbar; axis image; title('Original image');
% figure;imagesc(real(log10((im)))); colorbar; axis image; title('Original image');caxis([-4 0]);

soln_PD = asol;%/max(asol(:));
% soln = soln + soln_;
% figure; imagesc(log10(max(soln_, 1e-8))); colorbar, axis image; title('Recovered image');
% figure; imagesc(real(log10((soln_PD))));set(gca,'YDir','normal');axis off;
% colorbar, axis image; caxis([f_log_min f_log_max]);%title('Recovered image');
% colormap(cubehelix);


dirty = At(Gw' * cell2mat(ay));
dirty = 2 * real(dirty);
% normalized dirty image
dirtyn_ = dirty - min(dirty(:));
dirtyn_ = dirtyn_/max(dirtyn_(:));
dirtyn = dirtyn + dirtyn_;
% figure; imagesc(log10(max(dirtyn_, 1e-8))); colorbar, axis image; title('Dirty image');


% figure;imagesc(residualn), colorbar, axis image; title('Residual image');
% figure; imagesc(real(((residualn_PD))));set(gca,'YDir','normal');axis off;
% colorbar, axis image; caxis([-f_res f_res]);%title('Recovered image');
% colormap(cubehelix);

error_PD = im - soln_PD;
% figure; imagesc(real(((error_PD))));set(gca,'YDir','normal');axis off;
% colorbar, axis image; caxis([-f_err f_err]);%title('Recovered image');
% colormap(cubehelix);

SNR_PD = 20 * log10(norm(im(:))/norm(error_PD(:)));
TIME_PD=tend;

%% PPD
% run_pdfb_bpcon_par_sim_rescaled_precond

% PDFB parameter structure sent to the algorithm
param_pdfb_precond.im = im; % original image, used to compute the SNR
param_pdfb_precond.verbose = verbosity; % print log or not
param_pdfb_precond.nu1 = 1; % bound on the norm of the operator Psi
param_pdfb_precond.nu2 = evl_precond; % bound on the norm of the operator A*G
% param_pdfb_precond.gamma = 1e-6; % convergence parameter L1 (soft th parameter)
param_pdfb_precond.gamma = soft_threshold_values ; % convergence parameter L1 (soft th parameter)
% param_pdfb_precond.tau = 0.49; % forward descent step size
param_pdfb_precond.tau = forward_descent_step_size ; % forward descent step size
param_pdfb_precond.rel_obj = 0;1e-3; % stopping criterion
% param_pdfb_precond.max_iter = 20; % max number of iterations
param_pdfb_precond.max_iter = Maxiter ;
param_pdfb_precond.lambda0 = 1; % relaxation step for primal update
param_pdfb_precond.lambda1 = 1; % relaxation step for L1 dual update
param_pdfb_precond.lambda2 = 1; % relaxation step for L2 dual update
param_pdfb_precond.sol_steps = [inf]; % saves images at the given iterations

param_pdfb_precond.use_proj_elipse_fb = 1;
param_pdfb_precond.elipse_proj_max_iter = 200;
param_pdfb_precond.elipse_proj_min_iter = 1;
param_pdfb_precond.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid

param_pdfb_precond.use_reweight_steps = 0;
param_pdfb_precond.use_reweight_eps = 0;
param_pdfb_precond.reweight_steps = [500:100:2000 inf];
param_pdfb_precond.reweight_rel_obj = 1e-5; % criterion for performing reweighting
param_pdfb_precond.reweight_min_steps_rel_obj = 50;
param_pdfb_precond.reweight_alpha = 0.01; % omega^(0)
param_pdfb_precond.reweight_alpha_ff = 0.5; % exponential decay factor for omega^(0)


result_st = [];
result_st.sol = cell(num_tests, 1);
result_st.L1_v = cell(num_tests, 1);
result_st.L1_vp = cell(num_tests, 1);
result_st.L2_v = cell(num_tests, 1);
result_st.L2_vp = cell(num_tests, 1);
result_st.time = cell(num_tests, 1);
result_st.delta_v = cell(num_tests, 1);
result_st.sol_v = cell(num_tests, 1);
result_st.snr_v = cell(num_tests, 1);
result_st.sol_reweight_v = cell(num_tests, 1);
result_st.no_sub_itr_v = cell(num_tests, 1);
result_st.xcorr_v = cell(num_tests, 1);

result_st.snr = cell(num_tests, 1);
result_st.sparsity = cell(num_tests, 1);
result_st.no_itr = cell(num_tests, 1);


for i = 1:num_tests
    % wavelet mode is a global variable which does not get transfered
    % to the workes; we need to set it manually for each worker
    dwtmode('per');
    
    fprintf('Test run %i:\n', i);
    
    tstart_a = tic;
    fprintf(' Running pdfb_bpcon_par_sim_rescaled_precond\n');
    [result_st.sol{i}, result_st.L1_v{i}, result_st.L1_vp{i}, result_st.L2_v{i}, ...
        result_st.L2_vp{i}, result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}, result_st.no_sub_itr_v{i}, ~, ~, ...
        result_st.sol_reweight_v{i}] ...
        = pdfb_bpcon_par_sim_rescaled_precond(yT{i}, epsilonT{i}, epsilonTs{i}, epsilon{i}, epsilons{i}, A, At, T, aW, W, Psi, Psit, Psiw, Psitw, param_pdfb_precond);
    tend = toc(tstart_a);
    fprintf(' pdfb_bpcon_par_sim_rescaled_precond runtime: %ds\n\n', ceil(tend));
    
    result_st.time{i} = tend;
    if ~use_real_visibilities
        error = im - result_st.sol{i};
        result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
    end
    result_st.no_itr{i} = length(result_st.L1_v{i});
    
    wcoef = [];
    for q = 1:length(Psit)
        wcoef = [wcoef; Psit{q}(result_st.sol{i})];
    end
    result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
end

% results_prefix = 'pdfb_bpcon_par_sim_rescaled_precond';
% param_structure_name = 'param_pdfb_precond';

snr = 0;
sp = 0;
soln = zeros(Ny, Nx);
residualn = zeros(Ny, Nx);
dirtyn = zeros(Ny, Nx);
asol = result_st.sol{i};
ay = y{i};
aL1_v = result_st.L1_v{i};
aL2_v = result_st.L2_v{i};
aL1_vp = result_st.L1_vp{i};
aL2_vp = result_st.L2_vp{i};
% ay0 = y0{i};
aepsilon = epsilon{i};
aepsilonT = epsilonT{i};
R = length(ay);
ys = A(asol);

% results
% script_gen_figures;


% save data
% script_save_result_data;



asol = result_st.sol{i};
residual = At(Gw' * (cell2mat(ay) - Gw * ys));
residualn_PPD = real(residual)/evl;
% residualn = residualn + residualn_;

% figure;imagesc(log10(max(im, 1e-8))); colorbar; axis image; title('Original image');

soln_PPD = asol;%/max(asol(:));
% soln = soln + soln_;
% figure; imagesc(log10(max(soln_, 1e-8))); colorbar, axis image; title('Recovered image');
% figure; imagesc(log10(max(soln_PPD, 1e-4))); colorbar, axis image; title('Recovered image');caxis([-3.5 0]);
% colormap(cubehelix);
figure; imagesc(real(log10((soln_PPD))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([f_log_min f_log_max]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

% normalized dirty image
dirtyn_ = dirty - min(dirty(:));
dirtyn_ = dirtyn_/max(dirtyn_(:));
dirtyn = dirtyn + dirtyn_;
% figure; imagesc(log10(max(dirtyn_, 1e-8))); colorbar, axis image; title('Dirty image');

% figure;imagesc(residualn), colorbar, axis image; title('Residual image');
figure; imagesc(real(((residualn_PPD))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_res f_res]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

error_PDD = im - soln_PPD;
figure; imagesc(real(((error_PDD))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_err f_err]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

SNR_PPD = 20 * log10(norm(im(:))/norm(error_PDD(:)));
TIME_PPD=tend;

% %% ADMM
% % run_admm_bpconpar
% % ADMM parameter structure sent to the algorithm
% param_admm.im = im; % original image, used to compute the SNR
% param_admm.verbose = verbosity; % print log or not
% % param_admm.gamma = 1e-6*evl; % convergence parametter
% param_admm.gamma = soft_threshold_values_admm*evl; % convergence parametter
% param_admm.nu = evl; % bound on the norm of the operator A
% param_admm.rel_obj = 0;1e-4; % stopping criterion
% % param_admm.max_iter = 20; % max number of iterations
% param_admm.max_iter = Maxiter ;
% param_admm.tight_L1 = 0; % indicate if Psit is a tight frame (1) or not (0)
% param_admm.max_iter_L1 = 100; % max number of iterations to be performed for estimating the L1 prox
% param_admm.rel_obj_L1 = 1e-3; % stopping criterion for the L1 prox
% param_admm.pos_L1 = 1; % constrain for the positivity
% param_admm.nu_L1 = 1; % bound on the norm of the operator Psi
% param_admm.verbose_L1 = 0;  % print log or not
% param_admm.sol_steps = inf; % saves images at the given iterations
% 
% result_st = [];
% result_st.sol = cell(num_tests, 1);
% result_st.L1_v = cell(num_tests, 1);
% result_st.L1_vp = cell(num_tests, 1);
% result_st.L2_v = cell(num_tests, 1);
% result_st.L2_vp = cell(num_tests, 1);
% result_st.time = cell(num_tests, 1);
% result_st.delta_v = cell(num_tests, 1);
% result_st.sol_v = cell(num_tests, 1);
% result_st.snr_v = cell(num_tests, 1);
% 
% result_st.snr = cell(num_tests, 1);
% result_st.sparsity = cell(num_tests, 1);
% result_st.no_itr = cell(num_tests, 1);
% 
% for i = 1:num_tests
%     % wavelet mode is a global variable which does not get transfered
%     % to the workes; we need to set it manually for each worker
%     dwtmode('per');
%     
%     fprintf('Test run %i:\n', i);
%     
%     tstart_a = tic;
%     fprintf(' Running admm_bpconpar\n');
%     [result_st.sol{i}, ~, result_st.L1_v{i}, result_st.L2_v{i}, ...
%         result_st.delta_v{i}, result_st.sol_v{i}, result_st.snr_v{i}] ...
%         = admm_bpconpar(yT{i}, cell2mat(epsilonT{i}), cell2mat(epsilonTs{i}), A, At, T, W, Psiw, Psitw, param_admm);
%     tend = toc(tstart_a);
%     fprintf(' admm_bpconpar runtime: %ds\n\n', ceil(tend));
%     
%     result_st.time{i} = tend;
%     if ~use_real_visibilities
%         error = im - result_st.sol{i};
%         result_st.snr{i} = 20 * log10(norm(im(:))/norm(error(:)));
%     end
%     result_st.no_itr{i} = length(result_st.L1_v{i});
%     
%     wcoef = [];
%     for q = 1:length(Psit)
%         wcoef = [wcoef; Psit{q}(result_st.sol{i})];
%     end
%     result_st.sparsity{i} = sum(abs(wcoef) > 1e-3)/length(wcoef);
% end
% % results_prefix = 'admm_bpconpar';
% % param_structure_name = 'param_admm';
% 
% snr = 0;
% sp = 0;
% soln = zeros(Ny, Nx);
% residualn = zeros(Ny, Nx);
% dirtyn = zeros(Ny, Nx);
% asol = result_st.sol{i};
% ay = y{i};
% aL1_v = result_st.L1_v{i};
% aL2_v = result_st.L2_v{i};
% aL1_vp = result_st.L1_vp{i};
% aL2_vp = result_st.L2_vp{i};
% % ay0 = y0{i};
% aepsilon = epsilon{i};
% aepsilonT = epsilonT{i};
% R = length(ay);
% ys = A(asol);
% 
% % results
% % script_gen_figures;
% 
% % save data
% % script_save_result_data;
% 
% asol = result_st.sol{i};
% residual = At(Gw' * (cell2mat(ay) - Gw * ys));
% residualn_ADMM = real(residual)/evl;
% % residualn = residualn + residualn_;
% 
% % figure;imagesc(log10(max(im, 1e-8))); colorbar; axis image; title('Original image');
% 
% soln_ADMM = asol;%/max(asol(:));
% % soln = soln + soln_;
% % figure; imagesc(log10(max(soln_ADMM, 1e-8))); colorbar, axis image; title('Recovered image');caxis([-3.5 0]);
% % colormap(cubehelix);
% figure; imagesc(real(log10((soln_ADMM))));set(gca,'YDir','normal');axis off;
% colorbar, axis image; caxis([f_log_min f_log_max]);%title('Recovered image');
% colormap(cubehelix);
% set(gca,'position',[-0.05,0.02,0.98,0.934])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
% 
% % normalized dirty image
% dirtyn_ = dirty - min(dirty(:));
% dirtyn_ = dirtyn_/max(dirtyn_(:));
% dirtyn = dirtyn + dirtyn_;
% % figure; imagesc(log10(max(dirtyn_, 1e-8))); colorbar, axis image; title('Dirty image');
% 
% % figure;imagesc(residualn), colorbar, axis image; title('Residual image');
% figure; imagesc(real(((residualn_ADMM))));set(gca,'YDir','normal');axis off;
% colorbar, axis image; caxis([-f_res f_res]);%title('Recovered image');
% colormap(cubehelix);
% set(gca,'position',[-0.05,0.02,0.98,0.934])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
% 
% error_ADMM = im - soln_ADMM;
% figure; imagesc(real(((error_ADMM))));set(gca,'YDir','normal');axis off;
% colorbar, axis image; caxis([-f_err f_err]);%title('Recovered image');
% colormap(cubehelix);
% set(gca,'position',[-0.05,0.02,0.98,0.934])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
% 
% SNR_ADMM = 20 * log10(norm(im(:))/norm(error_ADMM(:)));
% TIME_ADMM =tend;

%% IUWT算法

addpath iuwt/


% dirty = At(Gw' * cell2mat(ay));
% dirty = 2 * real(dirty);
% % normalized dirty image
% dirtyn_ = dirty - min(dirty(:));
% dirtyn_ = dirtyn_/max(dirtyn_(:));
% dirtyn = dirtyn + dirtyn_;

y22=y{1,1}{1,1};
xishu=op_norm(@(x) Tw * A(x), @(x) At(Tw' * x), [Ny, Nx], 1e-4, 200, 1);
Dirtymap=dirtyn./xishu;%脏图


N= 1 ;%

positiveflg=1;
waveletflg=1;
level=4;

% 记录迭代过程中的ssim和rmse
times=0;
SSIM_total=[];
RMSE_total=[];
PSNR_total=[];


% [m,n]=find(PSF_iuwt==max(max(PSF_iuwt)));
% center=[m n];

    W_IUWT = @(x) IUWT(x,level);% W means IUWT  小波变换
    WT_IUWT = @(x) IUWT(x,-level);% WT  means inverse IUWT

[m,n]=size(Dirtymap);

%Computng the UV mask with the psf
% UV=fft2(circshift(PSF_iuwt,1-center));
% maxUV=max(max(abs(UV)));
% indx=find(abs(UV/maxUV)>0.05);
% invuv=zeros(m,n);
% invuv(indx)=1./(UV(indx).^2);

%Initialization
im_temp=Dirtymap;% 初始图像为脏图
X_temp=W_IUWT(zeros(m,n));% 小波变换域下的系数
weight=ones(size(X_temp));  %权重


lambda=lambda_IUWT;
gamma = gamma_IUWT ;% 正则化参数γ
iter_max  =  iter_max_IUWT ; 

tol_error = 1e-6;%停止误差
im_diff_err=1;%图像当前误差，初始为1
iter=1;%当前迭代次数
im_rec_hat=Dirtymap;
im_rec_last=Dirtymap;
t_new=1;
X_old=X_temp;

IUWT_SNR_huitu=[ ];
tstart_a = tic;
while (im_diff_err>=tol_error)&&(iter<=iter_max)  
    t_old=t_new;
%     im_rec=im_rec_hat + gamma*(At2(y{1,1}{1,1} - A2(im_rec_hat)))./eval;%Fu为U*F  采样矩阵乘以傅里叶变换矩阵
    
        temp1=A(im_rec_hat);
    temp2 = temp1(W{1,1});
    temp3=T{1,1}*temp2;
    im_rec=im_rec_hat + gamma*(At(Tw' *(y22 - temp3)))./xishu;%Fu为U*F  
    

    coefs=W_IUWT(im_rec);
            % Soft thresholding
        Shrink=abs(coefs)-lambda*weight;
        X_temp=sign(coefs).*((Shrink>0).*Shrink); % X_temp为τ收缩算子   论文中的公式（23）
        
        im_rec = WT_IUWT(X_temp); % thr_coefs为取阈值后的α

        %Updating t and X
        t_new=(1+sqrt(1+4*t_old^2))/2; % t_new为t（k+1）   论文中的公式（24）
        im_rec_hat=im_rec+(t_old-1)/t_new*(im_rec-im_rec_last);  % X为β（k+1），更新IUWT域中脏图的小波系数 维度为500*3500     论文中的公式（25）

%     im_rec = WT(thr_coefs); % thr_coefs为取阈值后的α
    
%     t_now = (1 + sqrt(1 + 4*t_last*t_last))/2;% t(k+1)
%     im_rec_hat = im_rec + (t_last - 1)/t_now*(im_rec - im_rec_last);%更新后的图像
     if iter>1
            im_diff_err = norm(im_rec_last(:) - im_rec(:),2)/norm(im_rec_last(:),2);%误差大小
     end
     % 更新im_rec x(k)、t_now t(k)
    im_rec_last = im_rec;
    t_old = t_new;
    
%     weight=1./(abs(W(im_rec))+.001);  % 权重
    weight=1./(abs(W_IUWT(im_rec))+.1);  % 权重
    iter = iter +1;
    
    temp=real(im_rec);
    %temp(temp<0)=0;

end
tend = toc(tstart_a); 

% im_rec(im_rec<0)=0;
u_IUWT=real(im_rec);

snr = 0;
sp = 0;
soln = zeros(Ny, Nx);
residualn = zeros(Ny, Nx);
dirtyn = zeros(Ny, Nx);

soln_IUWT = u_IUWT;%/max(asol(:));
% soln = soln + soln_;
% figure; imagesc(log10(max(soln_, 1e-8))); colorbar, axis image; title('Recovered image');
% figure;imagesc(real(log10((soln_)))); colorbar; axis image;title('Recovered image');caxis([-4 0]);
% figure;imagesc(real(log10((soln_WTF)))); colorbar; axis image;title('Recovered image');caxis([-3.5 0]);
% colormap(cubehelix);
figure; imagesc(real(log10((soln_IUWT))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([f_log_min f_log_max]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

% normalized dirty image
dirtyn_ = dirty - min(dirty(:));
dirtyn_ = dirtyn_/max(dirtyn_(:));
dirtyn = dirtyn + dirtyn_;
% figure; imagesc(log10(max(dirtyn_, 1e-8))); colorbar, axis image; title('Dirty image');

ys = A(soln_IUWT);
residual = At(Gw' * (cell2mat(ay) - Gw * ys));
residualn_IUWT = real(residual)/evl;
% residualn = residualn + residualn_;
% figure;imagesc(residualn), colorbar, axis image; title('Residual image');
figure; imagesc(real(((residualn_IUWT))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_res f_res]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

error_IUWT = im - u_IUWT;
figure; imagesc(real(((error_IUWT))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_err f_err]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

SNR_IUWT = 20 * log10(norm(im(:))/norm(error_IUWT(:)));
TIME_IUWT =tend;

%% 紧框架

% addpath('Others');
% addpath('Utilities');
% addpath('Toolbox\DFrT_Tensor');
% addpath('Solvers');
% 
% % evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-6, 200, verbosity);
% xishu=op_norm(@(x) Tw * A(x), @(x) At(Tw' * x), [Ny, Nx], 1e-4, 200, 1);
% 
% 
% % 使用极化SARA生成的加噪的可见度函数
% y22=y{1,1}{1,1};
% Dirtymap=dirtyn./xishu;%脏图
% 
% 
% % undersampling k-space
% % mask=real(UV);% mask采样矩阵
% 
% % row=1024;
% % column=1024;
% [row,column]=size(im);
% 
% % Fu = Fu_downsample(mask,row,column);
% 
% % y = Fu*Original;
% 
% 
% 
% % 初始化   
% % gamma = 0.75 ;
% gamma = gamma_WTF ;
% % iter_max  =  50 ; %最大迭代次数，传入参数niter
% iter_max  =  iter_max_WTF ;
% tol_error = 1e-6;%停止误差
% im_diff_err=1;%图像当前误差，初始为1
% iter=1;%当前迭代次数
% % img=Original;%原始图像，只用来计算ssim等，传入参数Original
% %Fu; % Fu为U*F  采样矩阵乘以傅里叶变换矩阵
% % imgMasked= Fu'*y;% mask掩模  Fu：k空间采样
% imgMasked=Dirtymap;
% 
% 
% nLvl  = 4 ; 
% % psi = TPCTFs(nLvl,6,row,column);
% psi = TPCTFs(nLvl,6,row,column);% 3,4,6
% im_rec_last = imgMasked;
% im_rec = imgMasked;
% im_rec_hat = im_rec;
% t_last = 1;
% 
% % temp1=A(im_rec_hat);
% % temp2 = temp1(W{1,1});
% % temp3=T{1,1}*temp2;
% 
% tstart_a = tic;
% TPCTF_SNR_huitu=[];
% while (im_diff_err>=tol_error)&&(iter<=iter_max)  
% %     im_rec=im_rec_hat + gamma*(At2(y22 - A2(im_rec_hat)))./xishu;%Fu为U*F  采样矩阵乘以傅里叶变换矩阵
% %     im_rec=im_rec_hat + gamma*(At(Tw' *(y22 - temp3)));%Fu
% %     im_rec=im_rec_hat + gamma*(At(Tw' *(y22 - temp3)))./xishu;%Fu为U*F 
% 
%     temp1=A(im_rec_hat);
%     temp2 = temp1(W{1,1});
%     temp3=T{1,1}*temp2;
%     im_rec=im_rec_hat + gamma*(At(Tw' *(y22 - temp3)))./xishu;%Fu为U*F  
% 
%     
%     coefs = psi*im_rec;% coefs为小波系数α，psi为紧框架D
%     thr_coefs =  thr_bishrink(coefs);% BS双变量收缩
%     im_rec = psi'*thr_coefs; % thr_coefs为取阈值后的α
%     t_now = (1 + sqrt(1 + 4*t_last*t_last))/2;% t(k+1)
%     im_rec_hat = im_rec + (t_last - 1)/t_now*(im_rec - im_rec_last);%更新后的图像
%      if iter>1
%             im_diff_err = norm(im_rec_last(:) - im_rec(:),2)/norm(im_rec_last(:),2);%误差大小
%      end
%      % 更新im_rec x(k)、t_now t(k)
%     im_rec_last = im_rec;
%     t_last = t_now;
%     
% %     im_rec_cir=(circshift(im_rec,1-center));
% %     figure;
% %     pcolor(xi(1,:),eta(:,1),real(im_rec_cir));xlabel('\it\xi');ylabel('\it\eta');title(strcat('结果'));
% %     shading interp;axis xy;colorbar;colormap(jet);%caxis([0 5]);
% %     c=colorbar;c.Label.String = 'log T/K';c.FontName = 'Times';
%     
%     %评判指标
% %     out.rlne_iter(iter) = RLNE(img,im_rec);
% %     out.ssim_iter(iter) = ssim(abs(img),abs(im_rec));
% %     out.psnr_iter(iter) = psnr(img,im_rec);
% %     fprintf('TPCTF-pFISTA ---> ssim: %.6f \t RLNE: %.6f \t  PSNR: %.6f \t #iter %d \n',out.ssim_iter(iter),out.rlne_iter(iter),out.psnr_iter(iter),iter);
%     iter = iter +1;
%     
% %     figure;
% %     pcolor(real(im_rec));title(strcat('im_rec'));
% %     shading interp;axis xy;colorbar;colormap(jet);
%     
% TPCTF_SNR=10*log10((sum(sum(im.^2)))/(sum(sum((im_rec-im).^2))));
% TPCTF_SNR_huitu=[TPCTF_SNR_huitu TPCTF_SNR];
% 
% end
% tend = toc(tstart_a); 
% 
% % figure;plot(TPCTF_SNR_huitu);
% 
% 
% u_MHWTF=im_rec;
% % %信噪比
% MHWT_SNR=10*log10((sum(sum(im.^2)))/(sum(sum((u_MHWTF-im).^2))));
% 
% snr = 0;
% sp = 0;
% soln = zeros(Ny, Nx);
% residualn = zeros(Ny, Nx);
% dirtyn = zeros(Ny, Nx);
% 
% asol = u_MHWTF;
% % residual = At(Gw' * (cell2mat(ay) - Gw * ys));
% % residualn_ = real(residual)/evl;
% % residualn = residualn + residualn_;
% 
% % figure;imagesc(log10(max(im, 1e-8))); colorbar; axis image; title('Original image');
% % figure;imagesc(real(log10((im)))); colorbar; axis image; title('Original image');caxis([-4 0]);
% 
% soln_WTF = asol;%/max(asol(:));
% % soln = soln + soln_;
% % figure; imagesc(log10(max(soln_, 1e-8))); colorbar, axis image; title('Recovered image');
% % figure;imagesc(real(log10((soln_)))); colorbar; axis image;title('Recovered image');caxis([-4 0]);
% % figure;imagesc(real(log10((soln_WTF)))); colorbar; axis image;title('Recovered image');caxis([-3.5 0]);
% % colormap(cubehelix);
% figure; imagesc(real(log10((soln_WTF))));set(gca,'YDir','normal');axis off;
% colorbar, axis image; caxis([f_log_min f_log_max]);%title('Recovered image');
% colormap(cubehelix);
% set(gca,'position',[-0.05,0.02,0.98,0.934])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
% 
% % normalized dirty image
% dirtyn_ = dirty - min(dirty(:));
% dirtyn_ = dirtyn_/max(dirtyn_(:));
% dirtyn = dirtyn + dirtyn_;
% % figure; imagesc(log10(max(dirtyn_, 1e-8))); colorbar, axis image; title('Dirty image');
% 
% ys = A(soln_WTF);
% residual = At(Gw' * (cell2mat(ay) - Gw * ys));
% residualn_WTF = real(residual)/evl;
% % residualn = residualn + residualn_;
% % figure;imagesc(residualn), colorbar, axis image; title('Residual image');
% figure; imagesc(real(((residualn_WTF))));set(gca,'YDir','normal');axis off;
% colorbar, axis image; caxis([-f_res f_res]);%title('Recovered image');
% colormap(cubehelix);
% set(gca,'position',[-0.05,0.02,0.98,0.934])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
% 
% error_WTF = im - u_MHWTF;
% figure; imagesc(real(((error_WTF))));set(gca,'YDir','normal');axis off;
% colorbar, axis image; caxis([-f_err f_err]);%title('Recovered image');
% colormap(cubehelix);
% set(gca,'position',[-0.05,0.02,0.98,0.934])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
% 
% SNR_WTF = 20 * log10(norm(im(:))/norm(error_WTF(:)));
% TIME_WTF =tend;

%% WTF

addpath('Others');
addpath('Utilities');
addpath('Toolbox\DFrT_Tensor');
addpath('Solvers');

% evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-6, 200, verbosity);
xishu=op_norm(@(x) Tw * A(x), @(x) At(Tw' * x), [Ny, Nx], 1e-4, 200, 1);


% 使用极化SARA生成的加噪的可见度函数
y22=y{1,1}{1,1};
Dirtymap=dirtyn./xishu;%脏图

% undersampling k-space
% mask=real(UV);% mask采样矩阵

[row,column]=size(im);

% Fu = Fu_downsample(mask,row,column);

% y = Fu*Original;


% 初始化   
gamma = gamma_WTF ;
iter_max  =  iter_max_WTF ;
tol_error = 1e-5; %停止误差
im_diff_err=1;%图像当前误差，初始为1
iter=1;%当前迭代次数
% img=Original;%原始图像，只用来计算ssim等，传入参数Original
%Fu; % Fu为U*F  采样矩阵乘以傅里叶变换矩阵
% imgMasked= Fu'*y;% mask掩模  Fu：k空间采样
imgMasked=Dirtymap;


nLvl  = 4 ; 
% psi = TPCTFs(nLvl,6,row,column);
psi = TPCTFs(nLvl,6,row,column);% 3,4,6
im_rec_last = imgMasked;
im_rec = imgMasked;
im_rec_hat = im_rec;
t_last = 1;

tor=0; %重启的参数
cnt = 0;
flag = 1;
p = 1 ;%1/2
q = 1 ;%1/2
r=4;

% temp1=A(im_rec_hat);
% temp2 = temp1(W{1,1});
% temp3=T{1,1}*temp2;

tstart_a = tic;
TPCTF_SNR_huitu=[];
while (im_diff_err>=tol_error)&&(iter<=iter_max)  
%     im_rec=im_rec_hat + gamma*(At2(y22 - A2(im_rec_hat)))./xishu;%Fu为U*F  采样矩阵乘以傅里叶变换矩阵
%     im_rec=im_rec_hat + gamma*(At(Tw' *(y22 - temp3)));%Fu
%     im_rec=im_rec_hat + gamma*(At(Tw' *(y22 - temp3)))./xishu;%Fu为U*F 
    im_rec_hat_old=im_rec_hat;%重启的参数

    temp1=A(im_rec_hat);
    temp2 = temp1(W{1,1});
    temp3=T{1,1}*temp2;
    im_rec=im_rec_hat + gamma*(At(Tw' *(y22 - temp3)))./xishu;%Fu为U*F  

    
    coefs = psi*im_rec;% coefs为小波系数α，psi为紧框架D
    thr_coefs =  thr_bishrink(coefs);% BS双变量收缩
    im_rec = psi'*thr_coefs; % thr_coefs为取阈值后的α
%     t_now = (1 + sqrt(1 + 4*t_last*t_last))/2;% t(k+1)
    t_now = (p + sqrt(q + r*t_last*t_last))/2;
    im_rec_hat = im_rec + (t_last - 1)/t_now*(im_rec - im_rec_last);%更新后的图像
     if iter>1
            im_diff_err = norm(im_rec_last(:) - im_rec(:),2)/norm(im_rec_last(:),2);%误差大小
     end
     
     % 重启自适应策略
    vk = (im_rec_hat_old(:)-im_rec(:))'*(im_rec(:)-im_rec_last(:));
    if vk >= -1e-3     % tor  -1e-3 
        cnt = cnt + 1;
%         if cnt>=4 % increase the value here if the condition number is big
        if cnt>=1
            if flag
                a_half = (4+1*((t_last-1) /t_now)) /5 ;
                xi = a_half^(1/30);
                flag = 0;
            end
            r = r * xi;
            if r<3.99
                t_lim = ( 2*p + sqrt( r*p^2 + (4-r)*q ) ) / (4 - r);
                t_now = max(2 * t_lim, t_now);
            end
        else
            t_now = 1;% restart
        end
        im_rec_hat = im_rec;% restart
    end
      
     % 更新im_rec x(k)、t_now t(k)
    im_rec_last = im_rec;
    t_last = t_now;
    
    iter = iter +1;

end
tend = toc(tstart_a); % time
im_rec(im_rec<0)=0; % enforce positivity
u_rec=real(im_rec); % real

u_MHWTF=u_rec;
% %信噪比
MHWT_SNR=10*log10((sum(sum(im.^2)))/(sum(sum((u_MHWTF-im).^2))));

snr = 0;
sp = 0;
soln = zeros(Ny, Nx);
residualn = zeros(Ny, Nx);
dirtyn = zeros(Ny, Nx);

asol = u_MHWTF;
% residual = At(Gw' * (cell2mat(ay) - Gw * ys));
% residualn_ = real(residual)/evl;
% residualn = residualn + residualn_;

% figure;imagesc(log10(max(im, 1e-8))); colorbar; axis image; title('Original image');
% figure;imagesc(real(log10((im)))); colorbar; axis image; title('Original image');caxis([-4 0]);

soln_WTF = asol;%/max(asol(:));
% soln = soln + soln_;
% figure; imagesc(log10(max(soln_, 1e-8))); colorbar, axis image; title('Recovered image');
% figure;imagesc(real(log10((soln_)))); colorbar; axis image;title('Recovered image');caxis([-4 0]);
% figure;imagesc(real(log10((soln_WTF)))); colorbar; axis image;title('Recovered image');caxis([-3.5 0]);
% colormap(cubehelix);
figure; imagesc(real(log10((soln_WTF))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([f_log_min f_log_max]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);


% normalized dirty image
dirtyn_ = dirty - min(dirty(:));
dirtyn_ = dirtyn_/max(dirtyn_(:));
dirtyn = dirtyn + dirtyn_;
% figure; imagesc(log10(max(dirtyn_, 1e-8))); colorbar, axis image; title('Dirty image');

ys = A(soln_WTF);
residual = At(Gw' * (cell2mat(ay) - Gw * ys));
residualn_WTF = real(residual)/evl;
% residualn = residualn + residualn_;
% figure;imagesc(residualn), colorbar, axis image; title('Residual image');

figure; imagesc(real(((residualn_WTF))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_res f_res]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

error_WTF = im - u_MHWTF;
figure; imagesc(real(((error_WTF))));set(gca,'YDir','normal');axis off;
colorbar, axis image; caxis([-f_err f_err]);%title('Recovered image');
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);

SNR_WTF = 20 * log10(norm(im(:))/norm(error_WTF(:)));
TIME_WTF =tend;


