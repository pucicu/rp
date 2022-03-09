function y = embed(varargin)
%EMBED   Creates embedding vector using time delay embedding
%     Y=EMBED(X,M,T) creates the embedding vector Y from the time
%     series X using a time delay embedding with dimension M and
%     delay T. The resulting embedding vector has length N-T*(M-1),
%     where N is the length of the original time series.
%
%     Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         y = embed(x,2,17); % embed into 2 dimensions using delay 17
%         plot(y(:,1),y(:,2))
%
%     Reference:
%         Packard, N. H., Crutchfield, J. P., Farmer, J. D., 
%         Shaw, R. S. (1980). Geometry from a time series. 
%         Physical Review Letters 45, 712-716.

% Copyright (c) 2012-2022
% Potsdam Institute for Climate Impact Research
% Norbert Marwan
% http://www.pik-potsdam.de
%
% $Date: 2021/11/22 12:34:39 $
% $Revision: 5.2 $
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
tau = 1; m = 1;

%% get input arguments
if nargin > 2
    tau = varargin{3};
end
if nargin > 1
    m = varargin{2};
end    

x = varargin{1}; % input vector
N = size(x,1) - (m-1)*tau; % length of embedding vector
if N <= 1 % check if row vector (transform it to column vector)
   x = x(:);
   N = length(x) - (m-1)*tau; % length of embedding vector
end

d = size(x,2); % number of columns
if d > N % check if long enough related to the number of columns
   error('Number of columns should be smaller than the length of the time series.')
end

%% create embedding vector 
if size(x,2) == 1 % input vector is one-column time series
    y = buffer(x,N,N-tau,'nodelay');
else % input vector is multi-column time series
    y = zeros(N, d*m);
    for i = 1:d
       y(:, (i-1)*m+[1:m] ) = buffer(x(:,i),N,N-tau,'nodelay');
    end
end
