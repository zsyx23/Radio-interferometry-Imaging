% Description:
%   Computes sky from given visibilities.
%
% Usage:
% 
%   [sky] = Vis2Sky_DirectFT(u,v,L,M,vis,Bmat)
% 
% Input parameters:
%   u,v         -   u,v sample points
%
%   L,M         -   l,m grid on which sky is defined
%
%   vis         -   Visibilities
%
% Output parameters:
%   sky         -   Vector image representing sky image
%
% Notes:
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

function [sky] = Vis2Sky_DirectFT(u,v,L,M,vis)

    % number of sky pixels to process at a time for memory-limited systems
    N_SKY_PATCH_SIZE = 512;

    N_pix = numel(L);
    N_vis = numel(u);
    sky = zeros(N_pix,1);
    N_runs = ceil(N_pix/N_SKY_PATCH_SIZE);
    for irun=1:N_runs
        pix_range = (((irun-1)*N_SKY_PATCH_SIZE+1):(irun*N_SKY_PATCH_SIZE))';
        if (pix_range(end) > N_pix)
            pix_range = pix_range(pix_range <= N_pix);
        end
        % exponential factor, i.e. non-uniform Fourier inverse transform matrix
        tfm_inv = exp(1i*2*pi*(L(pix_range)*u(:)' + M(pix_range)*v(:)'));
        % add this sky patch visibility
        sky(pix_range) = tfm_inv*vis/N_vis;
    end

end
