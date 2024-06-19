% Description:
%   Simulates visibilities for given u,v coverage and sky model using given
%   direction-dependent and baseline-dependent correction functions.
% 
% Usage:
% 
%   [vis] = Sky2Vis_DirectFT_BmatModel(u,v,L,M,sky,dd_terms,bd_terms)
% 
% Input parameters:
%   u,v         -   u,v sample points
%
%   L,M         -   l,m grid on which sky is defined
%
%   sky         -   Vector image representing sky model
%
%   dd_terms    -   Direction-dependent expansion functions.

%   bd_terms    -   Baseline-dependent expansion functions.
%
% Output parameters:
%   vis         -   Simulated visibilities
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

function [vis] = Sky2Vis_DirectFT_BmatModel(u,v,L,M,sky,dd_terms,bd_terms)
    
    % number of sky pixels to process at a time for memory-limited systems
    N_SKY_PATCH_SIZE = 512;
    
    N_pix = numel(L);
    N_vis = numel(u);
    vis = zeros(N_vis,1);
    N_runs = ceil(N_pix/N_SKY_PATCH_SIZE);
    N_terms = size(dd_terms,2);
    for iterm=1:N_terms
        % apply this term direction-dependent correction sky-plane correction
        this_sky = dd_terms(:,iterm) .* sky(:);
        this_vis = zeros(N_vis,1);
        for irun=1:N_runs
            pix_range = (((irun-1)*N_SKY_PATCH_SIZE+1):(irun*N_SKY_PATCH_SIZE))';
            if (pix_range(end) > N_pix)
                pix_range = pix_range(pix_range <= N_pix);
            end
            % exponential factor, including w-term
            tfm = exp(-1i*2*pi*(u(:)*L(pix_range)' + v(:)*M(pix_range)'));
            % add this sky patch visibility
            this_vis = this_vis + tfm*this_sky(pix_range);
        end
        % applly this term baseline-dependent correction visibility-plane correction
        this_vis = (bd_terms(iterm,:).') .* this_vis;
        % and add to the total visibility
        vis = vis + this_vis;
end
