% Description:
%   Simple check to determine if code is being run in Octave or Matlab 
% Usage:
% 
%   [tf] = is_octave()
% 
% Input parameters:
% 
% Output parameters:
%   tf      -   True: running in Octave, False: running in Matlab
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

function r = is_octave ()
    
    persistent x;
    
    if (isempty (x))
        x = exist ('OCTAVE_VERSION', 'builtin') == 5;
    end
    
    r = x;
    
end
