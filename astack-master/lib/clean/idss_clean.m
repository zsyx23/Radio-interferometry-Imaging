% Description:
%   Implementation of Image Domain Source Substraction (IDSS) CLEAN algorithm.
% 
% Usage:
% 
%   [clean,res,clean_conv,res_vec] = idss_clean(L,M,res_in,beam_basis,psf,psf_shape)
%   [clean,res,clean_conv,res_vec] = idss_clean(L,M,res_in,beam_basis,psf,psf_shape,TOL_SIGMA)
%   [clean,res,clean_conv,res_vec] = idss_clean(L,M,res_in,beam_basis,psf,psf_shape,TOL_SIGMA,LOOP_GAIN)
%   [clean,res,clean_conv,res_vec] = idss_clean(L,M,res_in,beam_basis,psf,psf_shape,TOL_SIGMA,LOOP_GAIN,MAX_ITER)
%   [clean,res,clean_conv,res_vec] = idss_clean(L,M,res_in,beam_basis,psf,psf_shape,TOL_SIGMA,LOOP_GAIN,MAX_ITER,PRINT_TO)
% 
% Input parameters:
%   L,M         -   l,m sample points
%
%   res_in      -   Input residual image with same shape as L and M (typically 
%                   the dirty image after backward transforming visibilities).
%
%   beam_basis  -   Matrix in which each column stores a DD basis function
%
%   psf         -   Matrix in which each column stores the PSF associated with the visibility 
%                   plane weighting applied to each DD basis function.
%
%   psf_shape   -   Dimensions of L,M plane over which PSFs are defined. Ideally large enough
%                   so that convolution with point-source near edge of image gives accurate
%                   PSF at diametrically opposite image edge.
% 
%   TOL_SIGMA   -   (Optional) Relative size of peak above residual std dev to be identified as
%                   source.
%
%   LOOP_GAIN   -   (Optional) Under-estimate source by this factor.
%
%   MAX_ITER    -   (Optional) Maximum number of iterations before forced termination.
%
%   PRINT_TO    -   (Optional) File handle for standard output.
% 
% Output parameters:
%   clean       -   Image containing point-sources
%  
%   res_out     -   Output residual image
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

function [clean,res_out,res_vec,iter] = idss_clean(L,M,res_in,beam_basis,psf,psf_shape,TOL_SIGMA,LOOP_GAIN,MAX_ITER,PRINT_TO)

    if (nargin < 7)
        TOL_SIGMA = 3; % look for peaks above TOL_SIGMA*stddev of residual
    end
    if (nargin < 8)
        LOOP_GAIN = 0.1; % point-source intensity scaling
    end
    if (nargin < 9)
        MAX_ITER = 500; % maximum number of iterations through main loop
    end
    FLAG_RES_VECTOR = 0;
    FLAG_NAN_IN_RES = 0;

    N_BACKSPACE = 14;
    if (nargout > 2)
        FLAG_RES_VECTOR = 1;
        res_vec = nan(MAX_ITER,1);
        N_BACKSPACE = N_BACKSPACE + 29;
    end

    % shape of image
    LM_shape = size(L);

    % number of basis functions
    N_modes = size(beam_basis,2);

    % initialise output residual to input residual
    res_out = res_in;

    % initialise other outputs
    clean = zeros(LM_shape);
    clean_conv = zeros(LM_shape);

    fprintf(PRINT_TO,'IDSS-CLEAN iteration %5d of %5d',0,MAX_ITER);
    if (nargout > 3)
        fprintf(PRINT_TO,' (residual norm = %10.3e)',0);
    end
    if (is_octave())
        fflush(PRINT_TO);
    end
    for iter=1:MAX_ITER
        if (mod(iter,floor(MAX_ITER/100)) == 0)
            fprintf(PRINT_TO,'%s%5d of %5d',repmat(char(8),[1 N_BACKSPACE]),iter,MAX_ITER);
        end
        [res_max,idx_max] = max(res_out(:));
        % test if peak is above residual noise
        res_std_dev = std(res_out(:));
        if (FLAG_RES_VECTOR)
            res_vec(iter) = norm((res_out(:)),2);
            if (mod(iter,floor(MAX_ITER/100)) == 0)
                fprintf(PRINT_TO,' (residual norm = %10.3e)',res_vec(iter));
                if (is_octave())
                    fflush(PRINT_TO);
                end
            end
        end
        if (res_max <= TOL_SIGMA*res_std_dev)
            fprintf(PRINT_TO,' ...residual image peak is below tolerance, terminating.')
            break;
        end
        % build combined PSF based on source location and DD basis function weights
        %cbeam = zeros(psf_shape);
        %for imode=1:N_modes
        %    cbeam = cbeam + reshape((beam_basis(idx_max,imode)*(psf(:,imode))),psf_shape);
        %end
        cbeam = reshape(beam_basis(idx_max,:)*(psf.'),psf_shape);
        ccomp = (LOOP_GAIN*res_max);
        % add clean component
        clean(idx_max) = clean(idx_max) + ccomp;
        % subtract combined PSF from residual
        image_sub = zeros(LM_shape);
        image_sub(idx_max) = ccomp;
        image_sub = conv2(image_sub,cbeam,'same');
        res_out = res_out - image_sub;
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

end

