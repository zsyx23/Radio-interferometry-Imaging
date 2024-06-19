% Description:
% 
% Usage:
% 
%   [u,v,w,N_vis,comb_idx] = get_uv_coverage(x,y,z,wavelength)
%   [u,v,w,N_vis,comb_idx] = get_uv_coverage(x,y,z,wavelength,ant_idx)
% 
% Input parameters:
%   x,y,z           -   Positions of antennas in [m].
%
%   wavelength      -   Wavelength at the operating frequency
%
%   ant_idx         -   (Optional) Index range for antennas to use.
% 
% Output parameters:
%   u,v,w           -   u,v,w coordinates for visibilities sampled by array
%                       (normalised to wavelength)
% 
%   N_vis           -   Number of visibilities
%
%   comb_idx        -   N_vis x 2 matrix in which each row gives the pair
%                       of antennas used to form the corresponding
%                       baseline.
%
% Notes:
%   Visibility function is sampled symmetrically, so that u,v,w will, for
%   every u0,v0,w0 also contain -u0,-v0,-w0 (except for the origin 0,0,0).
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

function [u,v,w,N_vis,comb_idx] = get_uv_coverage(x,y,z,wavelength,varargin)

    ant_idx = [];
    if (nargin > 4)
        ant_idx = varargin{1};
    end

    if (~isempty(ant_idx))
        x = x(ant_idx);
        y = y(ant_idx);
        z = z(ant_idx);
    end
    N_ant = numel(x);

    comb_idx = combnk(1:N_ant,2);
    u = diff(x(comb_idx),[],2)/wavelength;
    u = [0; u; -u];
    v = diff(y(comb_idx),[],2)/wavelength;
    v = [0; v; -v];
    w = diff(z(comb_idx),[],2)/wavelength;
    w = [0; w; -w];
    N_vis = numel(u);
    
end
