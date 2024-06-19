% Description:
%   Creates a random point-source distribution on the sky
% 
% Usage:
% 
%   [sky] = sim_sky_brightness_point_sources(N_pix,N_ptsrc,M_flux)
% 
% Input parameters:
%   N_pix       -   Number of pixels in l,m frame
%
%   N_ptsrc     -   Number of point sources
%
%   M_flux      -   Maximum flux density
% 
% Output parameters:
%   sky         -   Image vector containing point-source dirac deltas
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

function [sky] = sim_sky_brightness_point_sources(N_pix,N_ptsrc,M_flux)

    randI = randperm(N_pix);
    randI = randI(1:N_ptsrc);
    sky = zeros(N_pix,1);
    z_ptsrc = M_flux*rand(numel(randI),1);
    sky(randI) = z_ptsrc;

end
