% This script demonstrates how the forward calculation (sky to visibilities) may
% be computed using A-stacking by trading accuracy for computational efficiency.
% The main output (the last figure created) is essentially a reproduction of
% Fig.4 in [1].
%
% Intermediate results useful to other demonstration scripts are also stored:
%   The B-matrix in eq (8a) in [1].
%   The Direction-Dependent and Baseline-dependent correction factors, based on 
%   the SVD of B, see eqs (38) and (39) in [1].
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

fprintf(PRINT_TO,'Building B-matrix\n')
if (is_octave())
    fflush(PRINT_TO);
end

% Build exact baseline-dependent direction-dependent gains matrix 
B_exact = build_antenna_complex_cross_power_matrix(Field,baseline_idx,L,M,'x');
% Functions used later on expect a function handle that returns the BD-DD gains
% matrix for a given set of l,m coordinates. Use here a dummy-function that simply
% extracts the requested elements from the alread defined B-matrix.
BForward_exact = @(xL,xM) get_submatrix_of_B(L,M,B_exact,xL,xM);

% store the B-matrix for future calculations
if (is_octave())
    save -hdf5 B_matrix.hdf5 B_exact
else
    save B_matrix.mat B_exact
end

fprintf(PRINT_TO,'Simulating exact visibilities\n')
if (is_octave())
    fflush(PRINT_TO);
end

% Simulate exact visibilities, essentially an implementation of eq (7) in [1]
vis_exact = Sky2Vis_DirectFT_BmatHad(u,v,L,M,sky,BForward_exact);

fprintf(PRINT_TO,'Building BD-DD gain model (SVD of B-matrix)\n')
if (is_octave())
    fflush(PRINT_TO);
end

% Build BD-DD gain model based on SVD of the full B-matrix
[f_mat,sB,vB] = svd(B_exact.','econ'); % eq (38) [1]
sB = diag(sB);
clear vB; % alpha-coeffs can be computed using vB,sB or uB,B_exact
a_mat = (f_mat)'*(B_exact.'); % eq (39) in [1]

% store the DD-expansion functions (uB) and BD-weights (alpha_mat) for future
% calculations
if (is_octave())
    save -hdf5 BD_DD_model.hdf5 f_mat a_mat sB
else
    save BD_DD_model.mat f_mat a_mat sB
end

% Plot SV spectrum
figure;
sB_normalised_log = log10(sB/max(sB));
plot(sB_normalised_log,'b-x','linewidth',2);
grid on;
xlabel('Singular mode');
ylabel('Normalised singular value');
title('Singular value spectrum for B-matrix');

fprintf(PRINT_TO,'Computing visibilities for average DD gain\n')
if (is_octave())
    fflush(PRINT_TO);
end

% Compute visibilities with average DD gain, i.e. baseline-INdependent. The average
% DD gain for all antennas is on the first row of B-exact, and is included by 
% element-wise multiplication of the sky with that row before transforming to the
% visibility domain.
vis_onebeam = Sky2Vis_DirectFT(u,v,L,M, (sky .* B_exact(1,:).') );

% Compute l2-norm in visibility error vector, normalised to true vector l2-norm
err_onebeam = norm(vis_exact - vis_onebeam)/norm(vis_exact);

fprintf(PRINT_TO,'Computing visibilities for models of various levels of accuracy\n')
if (is_octave())
    fflush(PRINT_TO);
end

% Compute visibilities using BD-DD model of various numbers of terms, and calculate
% error for each case. Visibility compution is basically an implementation of
% eqs (18) in [1]
vis_model = zeros(size(vis_exact));
N_terms = numel(sB_normalised_log);
err_model = zeros(N_terms,1);
fprintf(PRINT_TO,'\t\tNumber of used terms: %3d of %3d',0,N_terms)
if (is_octave())
    fflush(PRINT_TO);
end
for iterm = 1:N_terms
%     fprintf(PRINT_TO,'%s%3d of %3d',repmat(char(8),[1,10]),iterm,N_terms)
    if (is_octave())
        fflush(PRINT_TO);
    end
    % apply direction-dependent correction term in sky plane
    this_sky = f_mat(:,iterm) .* sky(:); % eq (18a)
    % transform modified sky to visibilities withouth any DD effects
    this_vis = Sky2Vis_DirectFT(u,v,L,M,this_sky); % eq (18b)
    % apply baseline-dependent correction term in visibility plane
    vis_model = vis_model + a_mat(iterm,:).' .* this_vis(:); % eq (18d)
    % compute error relative to exact
    err_model(iterm) = norm(vis_exact - vis_model) / norm(vis_exact);
end
fprintf(PRINT_TO,'...done.\n')
if (is_octave())
    fflush(PRINT_TO);
end

% Plot the error in the visibilities as function of the number of terms in the
% BD-DD gain model
figure;
hold on;
plot(log10(err_model),'b-x','linewidth',2,'displayname','BD-DD model');
plot(xlim,log10(err_onebeam)*[1 1],'r--','linewidth',2,'displayname','Avg DD gain');
grid on;
xlabel('Number of terms');
ylabel('Normalised error [log10]');
title('Error in visibilities as function of number of terms in BD-DD model. Errof for avg DD gain also shown.');
ylim([-5 0]);
legend show;
