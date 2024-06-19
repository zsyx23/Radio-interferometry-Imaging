% Description:
%   Build the exact baseline-dependent direction-dependent gain matrix.
% 
% Usage:
% 
%   [B] = build_antenna_complex_cross_power_matrix(Field,comb_idx,L,M)
%   [B] = build_antenna_complex_cross_power_matrix(Field,comb_idx,L,M,pol)
% 
% Input parameters:
%   Field           -   CAESAR format structure that represents antenna
%                       voltage patterns. Should be in rectangular polarisation.
%
%   comb_idx        -   N_vis x 2 matrix in which each row is the pair of
%                       antennas used to form the corresponding baseline.
%
%   L,M             -   l,m coordinates for which effective areas are
%                       computed.
%
%   pol             -   (optional) Request polarisation. One of 'x', 'y',
%                       or 'z'.
% 
% Output parameters:
%   B               -   Antenna effective area matrix. Different
%                       polarisations are along third dimension. Baselines are
%                       ordered as: [u,v=0; u,v>0; u,v<0]. For u,v=0 the
%                       average power pattern over all antennas is used.
% 
% Notes:
%   The input Field is assumed to contain open-circuit voltages that appear
%   at antenna terminals for incident plane-waves with 1 V/m electric
%   field.
%
%   Construction of B uses cubic interpolation of Field data. NaNs
%   resulting from interpolation are replaced with nearest neighbour.
%
%   Phase reference for each pattern should be at the assumed position of
%   the antenna within the interferometer array.
%
%   comb_idx correspond to one half of symmetrically distributed u,v points, and
%   excludes the origin. The output is given for 2*N_vis + 1 baselines, including
%   the origin.
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

function [B] = build_antenna_complex_cross_power_matrix(Field,comb_idx,L,M,pol)

    if (nargin < 5)        
        pol = 'a'; % all polarisations is default
        pol_vec = [1 2 3];
    end
    
    N_vis = size(comb_idx,1);
    N_pix = numel(L);
    % get theta,phi corresponding to l,m
    theta = asin(sqrt(L(:).^2+M(:).^2));
%     phi = wrap_to_2pi(atan2(M(:),L(:)));
    phi = wrapTo2Pi(atan2(M(:),L(:))); %yx
    
    % interpolate to get field on l,m grid
    [N_theta,N_phi,N_pol,N_el] = size(Field.E);
    % if Field.PHI extends over [0,2*pi-dphi], add extra column
    dphi = mean(diff(Field.PHI(1,:)));
    if ( (360 - max(max(Field.PHI))*180/pi) > (0.5*dphi) )
        THETA = Field.THETA(:,[1:end 1]);
        PHI = [Field.PHI repmat(2*pi,[N_theta 1])];
        E = Field.E(:,[1:end 1],:,:);
    else
        THETA = Field.THETA;
        PHI = Field.PHI;
        E = Field.E;
    end
    if (strcmpi(pol,'x') || strcmpi(pol,'y') || strcmpi(pol,'z'))
        N_pol = 1;
        switch(lower(pol))
            case 'x'
                pol_vec = 1;
            case 'y'
                pol_vec = 2;
            case 'z'
                pol_vec = 3;
        end
    end
    E_on_lm = zeros(N_el,N_pix,N_pol);
    for ipol=pol_vec
        for iel=1:N_el
            E_on_lm(iel,:,ipol) = interp2(THETA.',PHI.',E(:,:,ipol,iel).',theta,phi,'cubic');
            % check for NaN and replace with nearest neighbour
            nan_idx = find(isnan(E_on_lm(iel,:,ipol)));
            if (~isempty(nan_idx))
                theta_nan = theta(nan_idx);
                phi_nan = phi(nan_idx);
                E_on_lm(iel,nan_idx,ipol) = interp2(THETA.',PHI.',E(:,:,ipol,iel).',theta_nan,phi_nan,'nearest');
            end
        end
    end
    
    % build B-matrix for one half of baselines
    B = zeros(N_vis,N_pix,N_pol);
    for ipol=pol_vec
        B(:,:,ipol) = E_on_lm(comb_idx(:,1),:,ipol).*conj(E_on_lm(comb_idx(:,2),:,ipol)); 
    end
    % build diagonal terms (using arithmetic mean of all antenna patterns)
    B0 = zeros(N_el,N_pix,N_pol);
    for ipol=pol_vec
        B0(:,:,ipol) = abs(E_on_lm(:,:,ipol).^2);
    end
    % put B-matrix components together for all u,v, including origin
    B = [mean(B0); B; conj(B)];
    
end
