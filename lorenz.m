function dy = lorenz(t,y,p1,p2)
% LORENZ   ODE for the Lorenz system
%     Usage with ode45:
%        [t x] = ode45('lorenz',[0 200],rand(1,3));
%     Changing default parameter r
%         r = 210;
%        [t x] = ode45('lorenz',[0 200],rand(1,3),[],r);

% Copyright (c) 2016-2019
% Potsdam Institute for Climate Impact Research, Germany
% Institute of Geosciences, University of Potsdam, Germany
% Norbert Marwan, K. Hauke Kraemer
% http://www.pik-potsdam.de
%
% This program is free software: you can redistribute it and/or modify it under the terms of the
% GNU Affero General Public License as published by the Free Software Foundation, either
% version 3 of the License, or (at your option) any later version.
% This program is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
% details.
% You should have received a copy of the GNU Affero General Public License along with this
% program. If not, see <http://www.gnu.org/licenses/>.

% parameter r 
if nargin<4
   r=28;
else
   r=p2;
end

dy = zeros(3,1);

% ODE of Lorenz 63
dy(1) = 10 * (y(2) - y(1));
dy(2) = r*y(1) - y(2) - y(1)*y(3);
dy(3) = y(1)*y(2) - (8/3)*y(3);

