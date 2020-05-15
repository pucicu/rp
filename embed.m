function y = embed(varargin)
%EMBED   Creates embedding vector using time delay embedding
%     Y=EMBED(X,M,T) creates the embedding vector Y from the time
%     series X using a time delay embedding with dimension M and
%     delay T. The resulting embedding vector has length N-T*(M-1),
%     where N is the length of the original time series.
%
%     Reference:
%         Packard, N. H., Crutchfield, J. P., Farmer, J. D., 
%         Shaw, R. S. (1980). Geometry from a time series. 
%         Physical Review Letters 45, 712-716.
%
%     Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         y = embed(x,2,17); % embed into 2 dimensions using delay 17
%         plot(y(:,1),y(:,2))

% Copyright (c) 2012-2019
% Potsdam Institute for Climate Impact Research
% Norbert Marwan
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
 

%% check input and output arguments
narginchk(1,3)
nargoutchk(0,1)

%% set default values for embedding parameters
t = 1; m = 1;

%% get input arguments
if nargin > 2
    t = varargin{3};
end
if nargin > 1
    m = varargin{2};
end    

x = varargin{1}(:); % we expect a column vector
Nx = length(x); % length of time series
NX = Nx-t*(m-1); % length of embedding vector


%% create embeeding vector 
% performed using MATLAB's vectorization
% result is an index series representing the time delay and dimension
for i = 1:m;
    jx(1+NX*(i-1):NX+NX*(i-1)) = 1+t*(i-1):NX+t*(i-1);
end

%% the final embedding vector
y = reshape(x(jx),NX,m);
