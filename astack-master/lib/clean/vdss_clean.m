% Description:
%   Implementation of Visibility Domain Source Subtraction (VDSS) CLEAN algorithm.
% 
% Usage:
% 
%   [clean,res,clean_conv,res_vec] = vdss_clean(L,M,vis,Sky2Vis,Vis2Sky)
%   [clean,res,clean_conv,res_vec] = vdss_clean(L,M,vis,Sky2Vis,Vis2Sky,TOL_SIGMA)
%   [clean,res,clean_conv,res_vec] = vdss_clean(L,M,vis,Sky2Vis,Vis2Sky,TOL_SIGMA,LOOP_GAIN)
%   [clean,res,clean_conv,res_vec] = vdss_clean(L,M,vis,Sky2Vis,Vis2Sky,TOL_SIGMA,LOOP_GAIN,MAX_ITER)
%   [clean,res,clean_conv,res_vec] = vdss_clean(L,M,vis,Sky2Vis,Vis2Sky,TOL_SIGMA,LOOP_GAIN,MAX_ITER,PRINT_TO)
% 
% Input parameters:
%   L,M         -   l,m sample points
%
%   vis         -   Visibilities
%
%   Sky2Vis     -   Function that takes sky image and returns visibilities.
%                   Sky model is given as N_pix x 1 vector and visibilities
%                   should be returned as N_vis x 1 vector.
%
%   Vis2Sky     -   Function that takes visibilities and returns sky image.
%                   Visibilities are given as N_vis x 1 vector and sky 
%                   image is returned as N_pix x 1 vector.
%
%   TOL_SIGMA   -   (Optional) Relative size of peak above residual std dev to be
%                   identified as source.
%
%   LOOP_GAIN   -   (Optional) Under-estimate source by this factor.
%
%   MAX_ITER    -   (Optional) Maximum number of iterations before forced 
%                   termination.
%
%   PRINT_TO    -   (Optional) File handle for standard output.
% 
% Output parameters:
%   clean       -   Image vector containing point-sources
%  
%   res         -   Residual image vector
%
%   clean_conv  -   Clean image convolved with non-zero dirty beam
%
%   res_vec     -   Vector containing l2-norm of vector representation of
%                   residual per iteration.
%
%   iter        -   Number of iterations performed.
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

function [clean,res,clean_conv,res_vec,iter] = vdss_clean(L,M,vis,Sky2Vis,Vis2Sky,TOL_SIGMA,LOOP_GAIN,MAX_ITER,PRINT_TO)

    if (nargin < 6)
        TOL_SIGMA = 1; % look for peaks above TOL_SIGMA*stddev of residual
    end
    if (nargin < 7)
        LOOP_GAIN = 0.3; % point-source intensity scaling
    end
    if (nargin < 8)
        MAX_ITER = 500; % maximum number of iterations through main loop
    end
    if (nargin < 9)
        PRINT_TO = stdout;
    end
    
    CONV_BEAM_MIN = 0.1;
    FLAG_RES_VECTOR = 0;
    FLAG_NAN_IN_RES = 0;

    N_BACKSPACE = 14;
    if (nargout > 3)
        FLAG_RES_VECTOR = 1;
        res_vec = nan(MAX_ITER,1);
        N_BACKSPACE = N_BACKSPACE + 29;
    end

    % shapes of data structures
    N_vis = numel(vis);
    LM_shape = size(L);

    % initial residual is dirty image
    res = Vis2Sky(vis);
    if (sum(isnan(res(:))) > 0)
        FLAG_NAN_IN_RES = 1;
        res(isnan(res)) = 0;
    end
    % get 2D representation of residual
    res = reshape(res,LM_shape);
    clean = zeros(LM_shape);
    % construct dirty beam
    weights = ones([N_vis,1]);
    dbeam = Vis2Sky(weights);
    % get 2D dirty beam
    dbeam = reshape(dbeam,LM_shape);
    [srcL,srcM,srcI,] = deal(nan(MAX_ITER,1));
    N_src = 0;
    fprintf(PRINT_TO,'VDSS-CLEAN iteration %5d of %5d',0,MAX_ITER);
    if (nargout > 3)
        fprintf(PRINT_TO,' (residual norm = %10.3e)',0);
    end
    if (is_octave())
        fflush(PRINT_TO);
    end
    for iter=1:MAX_ITER
        fprintf(PRINT_TO,'%s%5d of %5d',repmat(char(8),[1 N_BACKSPACE]),iter,MAX_ITER);
        % 2D find for maximum pixel in residual
        [max_pix,I_max] = max(real(res)); % force real part of residual
        [max_pix,J_max] = max(max_pix);
        I_max = I_max(J_max);
        %max_pix = real(max_pix);
        % test if peak is above residual noise
        res_std_dev = std(res(:));
        if (FLAG_RES_VECTOR)
            res_vec(iter) = norm(res(:),2);
            fprintf(PRINT_TO,' (residual norm = %10.3e)',res_vec(iter));
        end
        if (max_pix <= TOL_SIGMA*res_std_dev)
            fprintf(PRINT_TO,' ...residual image peak is below tolerance, terminating.')
            break;
        end
        % add new source / update existing source
        if (ismember(L(I_max,J_max),srcL) && ismember(M(I_max,J_max),srcM))
            this_src_idx = find(L(I_max,J_max) == srcL & M(I_max,J_max) == srcM);
            srcI(this_src_idx) = srcI(this_src_idx) + LOOP_GAIN*max_pix;%/M_dbeam;
        else
            N_src = N_src + 1;
            srcL(N_src) = L(I_max,J_max);
            srcM(N_src) = M(I_max,J_max);
            srcI(N_src) = LOOP_GAIN*max_pix;%/M_dbeam;
        end
        % add clean component
        clean(I_max,J_max) = clean(I_max,J_max) + LOOP_GAIN*max_pix;%/M_dbeam;
        % forward-calculation
        vis_clean = Sky2Vis(clean(:));
        vis_res = vis - vis_clean;
        % update residual
        res = Vis2Sky(vis_res);
        if (sum(isnan(res(:))) > 0)
            FLAG_NAN_IN_RES = 1;
            res(isnan(res)) = 0;
        end
        % get 2D representation of residual
        res = reshape(res,LM_shape);
        if (iter == MAX_ITER)
            fprintf(PRINT_TO,'\nWarning: Maximum number of iterations reached. Possibly unfound sources remaining.');
            break;
        end
        if (is_octave())
            fflush(PRINT_TO);
        end
    end
    fprintf(PRINT_TO,'\n');
    if (is_octave())
        fflush(PRINT_TO);
    end
    % convolve clean components with equivalent dirty beam
    dbeam_conv = dbeam;
    dbeam_conv(dbeam_conv < CONV_BEAM_MIN) = 0; % just remove non-zero parts
    clean_conv = conv2(clean,dbeam_conv,'same');
    % return vector representation of each output parameter
    [clean,res,clean_conv] = deal(clean(:),res(:),clean_conv(:));
    if (FLAG_RES_VECTOR)
        res_vec = res_vec(~isnan(res_vec));
    end

    if (FLAG_NAN_IN_RES)
        fprintf(PRINT_TO,'warning: NaNs found in residual (set to zero). Possible divide by zero in antenna beam?');
    end

    % get index range for a single patch
    function pix_range = get_pix_range(irun,N_SKY_PATCH_SIZE,N_pix)
        pix_range = (((irun-1)*N_SKY_PATCH_SIZE+1):(irun*N_SKY_PATCH_SIZE))';
        if (pix_range(end) > N_pix)
            pix_range = pix_range(pix_range <= N_pix);
        end
    end

end

