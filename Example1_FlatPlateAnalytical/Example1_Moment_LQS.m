function [m] = Example1_Moment_LQS(x,~,Par)
% This function gives the prediction of the linear quasi steady model for a flat plate

% By Igor Kavrakov

%%%%%%%%% COPYRIGHT NOTICE %%%%%%%%% 
%  This file is part of AeroGPBuff.
%  AeroGPBuff  is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  AeroGPBuff  is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with AeroGPBuff .  If not, see <https://www.gnu.org/licenses/>.

% Copyright (c) Igor Kavrakov, Guido Morgenthal, Allan McRobie 2024

%%  Mean Linear Quasi Steady of a Flat plate - Moment

h = tan(x(:,1).*Par.x_max); % Vert disp (h_prime / B)
a = x(:,2+Par.Lag).*Par.x_max; % alpha
a_d = x(:,4+2*Par.Lag).*Par.x_max; % a_prime (Aerodynamic centre)
w = tan(x(:,5+2*Par.Lag).*Par.x_max); % vert gust (w / U)
m = pi/2.*(a+h+1/4*a_d+w);      % Unormalised moment

m = m./Par.y_max;               % Normalised moment

end
