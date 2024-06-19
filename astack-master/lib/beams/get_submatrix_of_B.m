% Description:
%   Used to build B on l,m from given B0 over l0,m0.
% 
% Usage:
% 
%   [B] = get_submatrix_of_B(L0,M0,B0,L,M)
% 
% Input parameters:
%
%   L0,M0       -   Available sample points for B-matrix.
%
%   B0          -   B-matrix defined over L0,M0
%   
%   L,M         -   Requested subset of sample points for B-matrix.
% 
% Output parameters:
%
%   B           -   BD-DD gain matrix for the given l,m coordinates
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

function [B] = get_submatrix_of_B(L0,M0,B0,L,M)

    [tf,idx] = ismember([L(:) M(:)],[L0(:) M0(:)],'rows');
    B = B0(:,idx(tf==1));

end
