clear all
close all
clc

t1=clock;

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
% solve_HS = 1 ; % wide-band solver (rwLRJS) applied on all the channels
% solve_1B = 0; % single-channel solver (SARA) applied for each channel separately
solve_HS = 1 ;
solve_1B = 0 ; 
load_temp_result = 0;

%% Generate or Load ground-truth
if generate_groundtruth
    file_name = './simulated_data/data/W28_256.fits';
    Nx=256;
    Ny=256;
    c = 60;
%     f = linspace(1,2,c); % frequency bandwidth from 1 to 2 GHz with 15 channels
    f = linspace(1,2,c);
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
%     unds = 4; % take 1/unds images
    unds = 60 ;
    [x0,X0,f,c] = Generate_undersampled_cube(x0,f,Ny,Nx,c,unds);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% importdata W28_256.fits ;
% x0=double(ans);
% X0=x0(:);

% importdata 30dor_256.fits ;
% importdata M31_256.fits ;
% x0=double(ans);
% x0=x0./max(max(real(x0)));
% X0=x0(:);

% importdata W28.fits ;
% im=double(ans);
% im=im/(max(max(im)));

importdata CYGCBEST-1024.fits ;
im=double(ans);
im=im/(max(max(im)));
temp=zeros(1024,1024);
temp(275:751,:)=im;
im=temp;

[Nx,Ny]=size(im);
im=imresize(im,1024./Nx);
x0=double(im);
X0=x0(:);

figure; imagesc(real(log10((x0))));
set(gca,'YDir','normal');
axis off;
colorbar, axis image; caxis([-3.5,0]);
colormap(cubehelix);
% set(gca,'position',[-0.05,0.02,0.98,0.934])
% set(gcf,'unit','normalized','position',[0.2,0.2,0.26,0.4]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ch = [1 : c]; % number of channels loaded (note that this can be one).
clear c;

%% Generate data

% percentage = 0.3 ; [0.3 0.05 0.1 0.02 0.01 0.2]; % sampling rate
% percentage =  [0.02 0.05 0.1 0.02 0.01 0.2]; % sampling rate
percentage= 0.4 ;
for k = 1 : length(percentage)
    
    %diary(['diary_HI_60_' num2str(percentage(k))]);
    
    
    if generate_data
%         Generate_data
        % config parameters
        N = Nx * Ny;
        ox = 2; % oversampling factors for nufft
        oy = 2; % oversampling factors for nufft
        Kx = 8; % number of neighbours for nufft
        Ky = 8; % number of neighbours for nufft
        
        % preconditioning parameters
        param_precond.N = N; % number of pixels in the image
        param_precond.Nox = ox*Nx; % number of pixels in the image
        param_precond.Noy = oy*Ny; % number of pixels in the image
        param_precond.gen_uniform_weight_matrix = 1; %set weighting type
        param_precond.uniform_weight_sub_pixels = 1;
        
        
        % block structure
        
        regenerate_block_structure = 1;
        
        param_block_structure.use_density_partitioning = 0;
        param_block_structure.density_partitioning_no = 1;
        
        param_block_structure.use_uniform_partitioning = 0;
        param_block_structure.uniform_partitioning_no = 4;
        
        param_block_structure.use_equal_partitioning = 1;
%         param_block_structure.equal_partitioning_no = 4;
        param_block_structure.equal_partitioning_no = 4 ;
        
        param_block_structure.use_manual_frequency_partitioning = 0;
        % sparam.fpartition = [pi]; % partition (symetrically) of the data to nodes (frequency ranges)
        % sparam.fpartition = [0, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
        % sparam.fpartition = [-0.25*pi, 0, 0.25*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
        % sparam.fpartition = [-64/256*pi, 0, 64/256*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
        param_block_structure.fpartition = [icdf('norm', 0.25, 0, pi/4), 0, icdf('norm', 0.75, 0, pi/4), pi]; % partition (symetrically) of the data to nodes (frequency ranges)
        % sparam.fpartition = [-0.3*pi, -0.15*pi, 0, 0.15*pi, 0.3*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
        % sparam.fpartition = [-0.35*pi, -0.25*pi, -0.15*pi, 0, 0.15*pi, 0.25*pi, 0.35*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
        
        param_block_structure.use_manual_partitioning = 0;
        
        % samplig pattern parameters
        % options 'gaussian', 'file', 'gaussian+large-holes', 'file+undersample', ''gaussian+missing-pixels'
        % sampling_pattern = 'file+undersample';
        sampling_pattern = 'gaussian';
%         sampling_pattern = 'gaussian+large-holes';
        %  sampling_pattern = 'gaussian+missing-pixels';
        % sampling_pattern = 'file';
        
        % sparam.file_name = '/Volumes/Data/MeasSets/meerkat2h.ar.uvw.dat'; % file name for uv coverage
        sparam.file_name = '/Volumes/Data/MeasSets/ska.2h.ar.uvw.dat'; % file name for uv coverage
        sparam.p = percentage(k); % number of measurements as proportion of number of pixels to recover
        sparam.hole_number = 8000; % number of holes to introduce for 'gaussian+large-holes'
        sparam.hole_prob = 0.05; % probability of single pixel hole for 'gaussian+missing-pixels'
        sparam.hole_size = pi/60; % size of the missing frequency data
        sparam.fpartition = [pi]; % partition (symetrically) of the data to nodes (frequency ranges)
        % sparam.fpartition = [0, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
        % sparam.fpartition = [-0.25*pi, 0, 0.25*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
        % sparam.fpartition = [-0.3*pi, -0.15*pi, 0, 0.15*pi, 0.3*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
        % sparam.fpartition = [-0.35*pi, -0.25*pi, -0.15*pi, 0, 0.15*pi, 0.25*pi, 0.35*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
        sparam.sigma = pi/3; % variance of the gaussion over continous frequency
        sparam.sigma_holes = pi/3; % variance of the gaussion for the holes
        
        % generate the sampling pattern
        sparam.N = N; % number of pixels in the image
        sparam.Nox = ox*Nx; % number of pixels in the image
        sparam.Noy = oy*Ny; % number of pixels in the image
      
        %
        if generate_uv
            [u, v] = util_gen_sampling_pattern(sampling_pattern, sparam);
%             if save_data
%                 save(['./simulated_data/data/uv_60_HI_final_' num2str(percentage(k)) '.mat'],'-v7.3', 'u', 'v');
%             end
        else
            load(['./simulated_data/data/uv_60_HI_final_' num2str(percentage(k)) '.mat']);
            disp('coverage loaded successfully ;)')
        end
        
        figure;plot(u{1,1},v{1,1},'.r');
        
        % remove all the visibilities outside the range [-pi, pi]
        r = sqrt(u{1}.^2 + v{1}.^2);
        size(r(r>pi))
        bmax = max(r);
        
        u1 = u{1};
        v1 = v{1};
        [mm] = find(r > pi);
        u1(mm) = 0;
        v1(mm) = 0;
        
        u1 = u1(u1 ~= 0);
        v1 = v1(v1 ~= 0);
        
        r = sqrt(u1.^2 + v1.^2);
        size(r(r>pi))
        bmax = max(r);
        
%         u1 = u1/2;
%         v1 = v1/2;%默认除以2



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% load ('test_ska_240s.uvw.mat');
load ('test_vla_10s.uvw.mat');
uw=uvw(:,1);vw=uvw(:,2);
MAX=max(max(real(uvw)));
uw=uw./MAX*pi;
vw=vw./MAX*pi;
figure;plot(uw,vw,'.r');axis equal;

u1=uw; v1=vw;
clearvars uw; clearvars vw;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


        
        %scatter(u1,v1,'r.')
        
        %
        for i = ch
            
            uw{i} = (f(i)/f(1)) * u1;
            vw{i} = (f(i)/f(1)) * v1;
            
            % compute uniform weights (sampling density) for the preconditioning
            [aWw] = util_gen_preconditioning_matrix(uw{i}, vw{i}, param_precond);
            
            % set the weighting matrix, these are the natural weights for real data
            nWw = ones(length(uw{i}), 1);
            
            % set the blocks structure
            [u, v, ~, uvidx, aW{i}, nW] = util_gen_block_structure(uw{i}, vw{i}, aWw, nWw, param_block_structure);
            
            % measurement operator initialization
            fprintf('Initializing the NUFFT operator\n\n');
            tstart = tic;
            
            [A, At, G{i}, W{i}, Gw{i}] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
        end
        
        % Save data
%         if save_data
%             save('./simulated_data/data/data.mat','-v7.3', 'G', 'W', 'aW');
%         end
        
%         if save_full_operator && exist('Gw','var')
%             save('./simulated_data/data/Gw.mat','-v7.3', 'Gw');
%         end
        
        % Free memory
%         if free_memory
%             clear u v u1 v1 uw vw aWw nW nWw;
%         end

    end % generate data

    
%% 
% input_snr = 40; % noise level (on the measurements)
input_snr = 25 ;
num_tests = 1;
    

%% Generate Measurements

% definition for the stoping criterion
param_l2_ball.type = 'sigma';
param_l2_ball.sigma_ball = 2;

%
q = 1;
for m = 1 : length(input_snr)
    for n = 1 : num_tests
        
        [y0, y, Nm, sigma_noise] = util_gen_measurements(x0, G, W, A, input_snr(m));
        
        [epsilon,epsilons] = util_gen_data_fidelity_bounds(y, Nm, param_l2_ball, sigma_noise);
        
        y0_t{q} = y0;
        y_t{q} = y;
        epsilont{q} = epsilon;
        epsilons_t{q} = epsilons;
        
        q = q + 1;
        
    end
end
    
%% Solver 
%   Solver_simulated_data

% Compute full measurement operator spectral norm
if compute_Anorm
    F = afclean( @(x) HS_forward_operator_precond(x, Gw, A, aW));
    Ft = afclean( @(y) HS_adjoint_operator_precond(y, Gw, At, aW, Ny, Nx));
    Anorm = pow_method_op(F, Ft, [Ny Nx length(ch)]);
%     if save_data
%         save('./simulated_data/data/Anorm.mat','-v7.3', 'Anorm');
%     end
elseif load_data
    load('./simulated_data/data/Anorm.mat');
end

% Compute measurement operator spectral norm for each channel individually
if compute_Anorm_ch
    for i = ch
        Anorm_ch(i) = pow_method_op(@(x) sqrt(cell2mat(aW{i})) .* (Gw{i}*A(x)), @(x) At(Gw{i}' * (sqrt(cell2mat(aW{i})) .* x)), [Ny Nx 1]);
    end
%     if save_data
%         save('./simulated_data/data/Anorm_ch.mat','-v7.3', 'Anorm_ch');
%     end
elseif load_data
    load('./simulated_data/data/Anorm_ch.mat');
end

if exist('Gw','var')
    clear Gw;
end
% Generate initial epsilons by performing imaging with NNLS on each data block separately
if generate_eps_nnls
    % param_nnls.im = im; % original image, used to compute the SNR
    param_nnls.verbose = 2; % print log or not
    param_nnls.rel_obj = 1e-5; % stopping criterion
    param_nnls.max_iter = 1000; % max number of iterations
    param_nnls.sol_steps = [inf]; % saves images at the given iterations
    param_nnls.beta = 1;
    % solve nnls per block
    for i = ch
        eps_b{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            % printf('solving for band %i\n\n',i)
            [~,eps_b{i}{j}] = fb_nnls_blocks(yb{i}{j}, A, At, G{i}{j}, W{i}{j}, param_nnls);
        end
    end
%     if save_data
%         save('./simulated_data/data/eps.mat','-v7.3', 'eps_b');
%     end
elseif load_data
    load('./simulated_data/data/eps.mat');
end

% sparsity operator definition
nlevel = 4; % wavelet level
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used
% wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8'}; % wavelet basis to be used

[Psi1, Psit1] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
P = length(Psi1);

for k = 1 : P
    f = '@(y) HS_forward_sparsity(y,Psi1{';
    f = sprintf('%s%i},Ny,Nx);', f,k);
    Psi{k} = eval(f);
    
    ft = '@(x) HS_adjoint_sparsity(x,Psit1{';
    ft = sprintf('%s%i},%i);', ft,k,1);
    Psit{k} = eval(ft);
end

%
% [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
% Psi = @(y) HS_forward_sparsity(y,Psiw,Ny,Nx);
% Psit = @(x) HS_adjoint_sparsity(x,Psitw,length(wlt_basis));

% Pnorm = pow_method_op(Psit,Psi,[Ny Nx c]);

%
% Anorm1 = pow_method_op(@(x) Gw{1}*A(x), @(x) At(Gw{1}' * x), [Ny Nx c]);

% HSI parameter structure sent to the  HSI algorithm
param_HSI.verbose = 2; % print log or not
param_HSI.nu0 = 1; % bound on the norm of the Identity operator
param_HSI.nu1 = 1; % bound on the norm of the operator Psi
param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G
param_HSI.gamma0 = 1;
param_HSI.gamma = 1e-2;  %convergence parameter L1 (soft th parameter)
param_HSI.rel_obj = 5e-4; % stopping criterion
param_HSI.max_iter = 2000; % max number of iterations
param_HSI.use_adapt_eps = 0; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
param_HSI.adapt_eps_start = 200; % minimum num of iter before stating adjustment
param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
param_HSI.adapt_eps_tol_out = 1.001; % tolerance outside the l2 ball
param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
% param_HSI.adapt_eps_steps = 100;
param_HSI.adapt_eps_rel_obj = 5e-4; % bound on the relative change of the solution
param_HSI.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data

param_HSI.reweight_alpha = 1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
param_HSI.reweight_alpha_ff = 0.5;
param_HSI.total_reweights = 5; % -1 if you don't want reweighting
param_HSI.reweight_abs_of_max = 1; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised

param_HSI.use_reweight_steps = 1; % reweighting by fixed steps
param_HSI.reweight_step_size = 300; % reweighting step size
param_HSI.reweight_steps = [5000: param_HSI.reweight_step_size :10000];
param_HSI.step_flag = 1;

param_HSI.use_reweight_eps = 0; % reweighting w.r.t the relative change of the solution
param_HSI.reweight_max_reweight_itr = param_HSI.max_iter - param_HSI.reweight_step_size;
param_HSI.reweight_rel_obj = 5e-4; % criterion for performing reweighting
param_HSI.reweight_min_steps_rel_obj = 300; % min num of iter between reweights

param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
% param_HSI.elipse_proj_max_iter = 20 ;
param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid

%
for q = 1 : length(input_snr) * num_tests %number of tests x number of InSNRs
    
% L21 + Nuclear
if solve_HS
    [xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,iterh,rel_fval,nuclear,l21,norm_res,res] = pdfb_LRJS_Adapt_blocks_rwNL21_par_precond_new_sim(y_t{q}, epsilons_t{q}, A, At, aW, G, W, Psi, Psit, param_HSI, X0);
%     if save_data
%         save('./simulated_data/data/result.mat','-v7.3','xsol', 'v0', 'v1', 'v2', 'g', 'weights0', 'weights1', 'proj', 't_block','reweight_alpha', 'epsilon', 'iterh', 'rel_fval', 'nuclear', 'l21', 'norm_res', 'res');
%     end
    
    c = size(xsol,3);
    sol = reshape(xsol(:),numel(xsol(:))/c,c);
    SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)))
    psnrh = zeros(c,1);
    for i = ch
        psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
    end
    SNR_average = mean(psnrh)
    
end

% rwL11 Adaptive Blocks
if solve_1B
    for iii=1:length(ch)
        Xsol(:,:,iii)=zeros([Ny,Nx]);
    end
    for i = ch
        i;
%         param_HSI.nu2 = Anorm_ch(i); % bound on the norm of the operator A*G
        param_HSI.nu2 = Anorm;
        [xsol,v1,v2,g,weights1,proj,t_block,reweight_alpha,epsilon,iterh,rel_fval,l11,norm_res,res] = pdfb_L11_Adapt_blocks_rw_par_precond_new({y_t{q}{i}},{epsilons_t{q}{i}}, A, At, {aW{i}}, {G{i}}, {W{i}}, Psi, Psit, param_HSI, i);
        Xsol(:,:,i)=xsol;
        
        c = size(Xsol,3);
        sol = reshape(Xsol(:),numel(Xsol(:))/c,c);
        SNR = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
        psnrh = zeros(c,1);
        for i = ch
            psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
        end
        SNR_average = mean(psnrh)
        
%         if save_data
%             save(['./simulated_data/data/result_b' num2str(i) '.mat'],'-v7.3', 'xsol', 'v1', 'v2', 'g', 'weights1', 'proj', 't_block','reweight_alpha', 'epsilon', 'iterh', 'rel_fval', 'l11', 'norm_res', 'res');
%         end
    end
    xsol=Xsol;
end

end
    
end


% figure;
% pcolor(real(log10(x0(:,:,1))));title('x0');
% shading interp;axis xy;colorbar;colormap(jet); caxis([-3.5,1]);

% 结果
figure; imagesc(real(log10((xsol(275:751,:)))));
set(gca,'YDir','normal');
axis off;
colorbar, axis image; caxis([-3.5,0]);
colormap(cubehelix);
set(gca,'position',[-0.05,0.02,0.98,0.934])
set(gcf,'unit','normalized','position',[0.2,0.2,0.4,0.3]);



