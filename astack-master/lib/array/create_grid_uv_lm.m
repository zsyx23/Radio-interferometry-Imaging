% Description:
%   Creates a sensibly sized l,m grid and corresponding FFT u,v grid based
%   on the given u,v sampling function.
% 
% Usage:
% 
%   [L,M,N_pix,U,V,N_grid_vis] = create_grid_uv_lm(u,v,M_lm,OVER_SAMPLE)
% 
% Input parameters:
%   u,v         -   u,v sample points
%
%   M_lm        -   Maximum extent of l,m frame
%
%   OVER_SAMPLE -   Over-sample factor as a power of 2
% 
% Output parameters:
%   L,M         -   l,m grid
%
%   N_pix       -   Number of pixels in l,m grid
%
%   U,V,        -   u,v on FFT grid
%
%   N_grid_vis  -   Number of visibilities in u,v grid
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

function [L,M,N_pix,U,V,N_grid_vis] = create_grid_uv_lm(u,v,M_lm,OVER_SAMPLE)

    M_uv = max(sqrt(abs(u).^2 + abs(v).^2)); % roughly extents of sampled visibility
    m_lm = 1/(M_uv*2); % get sensible resolution in lm-plane from extents of sampled visibility
    N_FFT = 2^(nextpow2(ceil(M_lm/m_lm)) + OVER_SAMPLE); % determine FFT-size
    m_lm = 2*M_lm/N_FFT; % update lm-plane resolution to make perfect N_FFT*m_lm fitt
    [l,m] = deal(-M_lm:m_lm:(M_lm-m_lm));
    [L,M] = meshgrid(l,m);
    M_uv = (1/m_lm)/2; % udpate uv-plane extents to match lm-resolution
    m_uv = 2*M_uv/N_FFT; % size of gridded uv-plane should be same as lm-plane
    [u_gridded,v_gridded] = deal(-M_uv:m_uv:(M_uv-m_uv));
    [U,V] = meshgrid(u_gridded,v_gridded);
    N_pix = numel(L);
    N_grid_vis = numel(U);        
        
end
