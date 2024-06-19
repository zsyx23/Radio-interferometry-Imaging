% Description:
%   Produce map from the given image vector.
% 
% Usage:
% 
%   [] = plot_sky_map(L,M,I)
% 
% Input parameters:
%   L,M         -   l,m coordinates for pixels
%
%   I           -   Image vector
% 
% Output parameters:
%   hf          -   Handle to figure object.
%
%   ha          -   Handle to axes object
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

function [hf,ha] = plot_sky_map(L,M,I,title_str)

    % check real/imaginary relation
    TOL_IMAG_OF_REAL =1e-9;
    imag_sum = sum(abs(imag(I(:))));
    real_sum = sum(abs(real(I(:))));
    if (imag_sum/real_sum > TOL_IMAG_OF_REAL)
        warning('Significant imaginary component detected in image vector.');
    end
    I_map = reshape(I,size(L));
    hf = figure;
    ha = axes;
    hold on;
%     set(pcolor2(L,M,real(I_map)),'LineStyle','none');
    set(pcolor(L,M,real(I_map)),'LineStyle','none'); %yx
    axis image;%(1.1*[-1 1 -1 1]*max([L(:); M(:)]));
    if (nargin < 4)
        title_str = 'Intensity map';
    end
    title(title_str);
    xlabel('l');
    ylabel('m');
    colorbar;

end
