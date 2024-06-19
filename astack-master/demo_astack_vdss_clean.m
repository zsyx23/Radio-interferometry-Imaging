% This script demonstrates the visibility domain source subtraction (VDSS) CLEAN
% algorithm using A-stacking. The main output (the last figure created) is 
% essentially a reproduction of Fig.5 in [1].
%
% This script loads data stored from demo_astack_forward_calculation.m:
%   B_matrix -- B-matrix (BD-DD gain)
%   f_mat, a_mat -- DD expansion functions and BD weights
% Note that these variables depend on the image plane and visibility plane
% discretisations.
%
% This script was developed in Octave (3.8.2), but may also run (possibly after
% minor tweaks) in Matlab. Under octave, the following octave-forge packages
% are required:
%   statistics
%
% References:
%   [1] A. Young, S.J. Wijnholds, T.D. Carozzi, R. Maaskant, M.V. Ivashina, and 
%   D.B. Davidson, ``Efficient correction for both direction-dependent and 
%   baseline-dependent effects in interferometric imaging: An A-stacking 
%   framework,'' Astronomy & Astrophysics, 2015.
%
% Copyright (C) 2014  Andre Young (ayoung.at.sun@gmail.com)
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
%   MA 02110-1301, USA.
%

% clear workspace
clear all
close all
clc

% path to libraries
addpath(genpath('./lib'));

% physical constants
c0 = 299792458;

PRINT_TO = stdout;

fprintf(PRINT_TO,'Loading station data\n')
if (is_octave())
    fflush(PRINT_TO);
end

% Load antenna data:
%   Field -- Structure that holds antenna voltage patterns. The structure has
%   the following fields:
%       .THETA -- 2D array that contains theta-coordinates of field samples (in
%       radians). Measured from theta=0 at zenith.
%       .PHI -- 2D array that contains phi-coordinates of field samples (in
%       radians). Measured from phi=0 along the x-axis.
%       .E -- 4D array that contains complex-valued voltages pattern samples.
%       The first two dimensions correspond to entries in .THETA and .PHI, the
%       third dimesion determines polarisation (x, y, and z components), and the
%       fourth dimension determines the antenna in the array. X-oriented dipoles
%       are entries 1:96.
%       .Polarization -- Indicates the polarisation convention, 'rectangular' for
%       Cartesian coordinates in this case.
%   freq -- The frequency at which patterns were obtained.
%   x,y,z -- Coordinates where antennas are located.
if (is_octave())
    load('./dat/LOFAR_LBA_AntennaData.hdf5');
else
    % Data was saved in Matlab v6 binary format from within Octave.
    load('./dat/LOFAR_LBA_AntennaData.mat');
end
Field.E = Field.E(:,:,:,1:96); % only use x-dipoles

% get u,v for snapshot towards zenith
wavelength = c0/freq;
[u,v,w,N_vis,baseline_idx] = get_uv_coverage(x,y,z,wavelength);
N_bl_oneside = size(baseline_idx,1);
N_bl_oneside_origin = N_bl_oneside+1;
N_bl_twoside = N_bl_oneside+N_bl_oneside_origin;
% simplify by setting w=0
w = zeros(size(w));

% setup l,m <--> u,v FFT grid
M_lm = 0.5; % image extents in lm-plane
N_FFT_OVER_SAMPLE = 0; % over sample factor (visibilities get zero-padded)
[L,M,N_pix,U_gridded,V_gridded,N_grid_vis] = create_grid_uv_lm(u,v,M_lm,N_FFT_OVER_SAMPLE);

fprintf(PRINT_TO,'Generating sky model\n')
if (is_octave())
    fflush(PRINT_TO);
end

% Generate a toy sky model over the l,m grid. The sky consists of a number of 
% point-like sources and otherwise empty sky.
N_ptsrc = 7; % number of point-like source
M_flux = 100; % maximum flux in arbitrary units
sky = sim_sky_brightness_point_sources(N_pix,N_ptsrc,M_flux);

% Visualise the sky model
plot_sky_map(L,M,sky,'True sky');

fprintf(PRINT_TO,'Loading B-matrix\n')
if (is_octave())
    fflush(PRINT_TO);
end
% Load B-matrix from file.
if (is_octave())
    load B_matrix.hdf5
else
    load B_matrix.mat
end
% Define function that returns B-matrix for given sky coordinates
BForward_exact = @(xL,xM) get_submatrix_of_B(L,M,B_exact,xL,xM);

% Simulate measured visibilities (exclude noise)
fprintf(PRINT_TO,'Simulate measured visibilities\n')
if (is_octave())
    fflush(PRINT_TO);
end
vis_exact = Sky2Vis_DirectFT_BmatHad(u,v,L,M,sky,BForward_exact);

% Load BD-DD gain model from file
fprintf(PRINT_TO,'Loading BD-DD gain model\n')
if (is_octave())
    fflush(PRINT_TO);
end
if (is_octave())
    load BD_DD_model.hdf5
else
    load BD_DD_model.mat
end
sB_normalised_log = log10(sB/max(sB));

% Set parameters for VDSS-CLEAN implementation
TOL_SIGMA = 1;
LOOP_GAIN = 0.5;
MAX_ITER = 100;

% For running VDSS-CLEAN we need to define functions that do the forward and 
% backward calculations. For the backward calculation, visibilities to sky, we 
% compute the dirty image and apply an element-wise direction-dependent correction 
% through dividing by the average antenna gain pattern.
Vis2Sky = @(x_vis) Vis2Sky_DirectFT(u,v,L,M,x_vis) ./ (B_exact(1,:).');

% For the forward calculation we define a different number of functions, depending
% on which model we use for DD-effects. This one uses the average DD gain, i.e.
% baseline-INdependent gain
Sky2Vis_onebeam = @(x_sky) Sky2Vis_DirectFT(u,v,L,M, (x_sky .* (B_exact(1,:).')) );
fprintf(PRINT_TO,'Average beam (baseline-independent) DD gain\n');
if (is_octave())
    fflush(stdout);
end
[clean_onebeam,res_onebeam,clean_conv_onebeam,res_vec_onebeam,iter_onebeam] = ...
    vdss_clean(L,M,vis_exact,Sky2Vis_onebeam,Vis2Sky,TOL_SIGMA,LOOP_GAIN,MAX_ITER,PRINT_TO);
    
% This one uses BD-DD gain model with all terms for which singular values are
% above 1/100 of maximum.
N_terms_2 = find(sB_normalised_log > -2,1,'last');
% The following line of code uses eq (18) in [1] directly to do the forward 
% calculation:
%   Sky2Vis_model_2 = @(x_sky) Sky2Vis_DirectFT_BmatModel(u,v,L,M,x_sky,f_mat(:,1:N_terms_2),a_mat(1:N_terms_2,:));
% but if sufficient storage is available in memory, then we can also limit the
% number of terms in the model by calculating the reduced-rank matrix as follows:
B_reduced_rank = (f_mat(:,1:N_terms_2) * a_mat(1:N_terms_2,:)).';
% and then using this reduced-rank matrix in eq (7) in [1]. For more realistic
% problems this would typically not be an option.
BForward_model = @(xL,xM) get_submatrix_of_B(L,M,B_reduced_rank,xL,xM);
Sky2Vis_model = @(x_sky) Sky2Vis_DirectFT_BmatHad(u,v,L,M,x_sky,BForward_model);
fprintf(PRINT_TO,'BD-DD gain model with %d terms (sigma > 0.01)\n',N_terms_2);
if (is_octave())
    fflush(stdout);
end
[clean_model_2,res_model_2,clean_conv_model_2,res_vec_model_2,iter_model_2] = ...
    vdss_clean(L,M,vis_exact,Sky2Vis_model,Vis2Sky,TOL_SIGMA,LOOP_GAIN,MAX_ITER,PRINT_TO);

% This one uses BD-DD gain model with all terms for which singular values are
% above 1/1000 of maximum.
N_terms_3 = find(sB_normalised_log > -3,1,'last');
%   Sky2Vis_model_3 = @(x_sky) Sky2Vis_DirectFT_BmatModel(u,v,L,M,x_sky,f_mat(:,1:N_terms_3),a_mat(1:N_terms_3,:));
B_reduced_rank = (f_mat(:,1:N_terms_3) * a_mat(1:N_terms_3,:)).';
BForward_model = @(xL,xM) get_submatrix_of_B(L,M,B_reduced_rank,xL,xM);
Sky2Vis_model = @(x_sky) Sky2Vis_DirectFT_BmatHad(u,v,L,M,x_sky,BForward_model);
fprintf(PRINT_TO,'BD-DD gain model with %d terms (sigma > 0.001)\n',N_terms_3);
if (is_octave())
    fflush(stdout);
end
[clean_model_3,res_model_3,clean_conv_model_3,res_vec_model_3,iter_model_3] = ...
    vdss_clean(L,M,vis_exact,Sky2Vis_model,Vis2Sky,TOL_SIGMA,LOOP_GAIN,MAX_ITER,PRINT_TO);

% This one uses BD-DD gain model with as many terms as antennas in the array
N_terms_ant = size(Field.E,4);
%   Sky2Vis_model_ant = @(x_sky) Sky2Vis_DirectFT_BmatModel(u,v,L,M,x_sky,f_mat(:,1:N_terms_ant),a_mat(1:N_terms_ant,:));
B_reduced_rank = (f_mat(:,1:N_terms_ant) * a_mat(1:N_terms_ant,:)).';
BForward_model = @(xL,xM) get_submatrix_of_B(L,M,B_reduced_rank,xL,xM);
Sky2Vis_model = @(x_sky) Sky2Vis_DirectFT_BmatHad(u,v,L,M,x_sky,BForward_model);
fprintf(PRINT_TO,'BD-DD gain model with %d terms (N_ant)\n',N_terms_ant);
if (is_octave())
    fflush(stdout);
end
[clean_model_ant,res_model_ant,clean_conv_model_ant,res_vec_model_ant,iter_model_ant] = ...
    vdss_clean(L,M,vis_exact,Sky2Vis_model,Vis2Sky,TOL_SIGMA,LOOP_GAIN,MAX_ITER,PRINT_TO);

% Finally, do VDSS-CLEAN using the exact BD-DD gain matrix
Sky2Vis_exact = @(x_sky) Sky2Vis_DirectFT_BmatHad(u,v,L,M,x_sky,BForward_exact);
fprintf(PRINT_TO,'Exact BD-DD gain\n');
if (is_octave())
    fflush(stdout);
end
[clean_exact,res_exact,clean_conv_exact,res_vec_exact,iter_exact] = ...
    vdss_clean(L,M,vis_exact,Sky2Vis_exact,Vis2Sky,TOL_SIGMA,LOOP_GAIN,MAX_ITER,PRINT_TO);

% Plot residual decrease over iterations for various solutions
figure;
axes;
hold on;
title(sprintf('VDSS-CLEAN convergence (tol = %d, loop gain = %05.3f, max iter = %d)',TOL_SIGMA,LOOP_GAIN,MAX_ITER));
xlabel('Number of iterations');
ylabel('Residual l2-norm');
plot(log10(res_vec_onebeam),'r','LineWidth',2,'DisplayName','Avg beam');
plot(log10(res_vec_model_2),'b','LineWidth',2,'DisplayName',sprintf('N_b = %d (sigma > 0.01)',N_terms_2));
plot(log10(res_vec_model_3),'g','LineWidth',2,'DisplayName',sprintf('N_b = %d (sigma > 0.001)',N_terms_3));
plot(log10(res_vec_model_ant),'k','LineWidth',2,'DisplayName','N_b = N_{ant} = 96');
plot(log10(res_vec_exact),'m','LineWidth',2,'DisplayName','Exact beam');
legend show
grid on
