%% MATLAB Example for calculating perpendicular and iso-directional RP

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


%% calculate example system Lorenz
[t x] = ode45('lorenz',[0 200],rand(1,3));
%x = resample(x(10300:11000,:),4,1);
x = x(10300:11000,:);

% perpendicular RP
[R1, SP, R0] = rp_perp(x,5,.5);

% iso-directional RP
R2 = rp_iso(x,5,.01);

subplot(131)
imagesc(R0)
title(sprintf('RR=%.3f',mean(R0(:))))

subplot(132)
imagesc(R1)
title(sprintf('RR=%.3f',mean(R1(:))))

subplot(133)
imagesc(R2)
title(sprintf('RR=%.3f',mean(R2(:))))


%% calculate example sine
xSin = sin(linspace(0,2*pi*7,1000));
x = embed(xSin,2,36);

% perpendicular RP (should be empty)
[R1, SP, R0] = rp_perp(x,.5,.1);

% iso-directional RP
R2 = rp_iso(x,.5,.001);


subplot(131)
imagesc(R0)

subplot(132)
imagesc(R1)

subplot(133)
imagesc(R2)
