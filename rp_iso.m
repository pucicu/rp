function [iso_r,sp,r] = rp_iso(varargin)
% RP_ISO    Calculates the isodirectional recurrence plot
%    R=RP_ISO(X,E,W) calculates the isodirectional recurrence plot R 
%    from an embedding vector X and using the threshold E for the
%    vector distances and threshold W for the angle to be 
%    considered as isodirectional.
%
%    R=RP_ISO(X,E,W,TAU) estimates tangential vector using time delay TAU.
%
%    Example:
%        [t x] = ode45('lorenz',[0 100],[-6.2 -10 14]);
%        [R1, SP, R0] = rp_iso(x(3000:5000,:),10,.2);
%
%        nexttile
%        imagesc(R0) % regular RP
%        axis square
%
%        nexttile
%        imagesc(R1) % isodirectional RP
%        axis square

% Copyright (c) 2016-2019
% Potsdam Institute for Climate Impact Research
% K. Hauke Kraemer, Norbert Marwan
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

%% check input
narginchk(1,5)
nargoutchk(0,3)

%% set default values for input parameters
e = 1; % recurrence threshold
w = .001; % angle threshold
tau = 1; % distance between points at phase space vector for constructing the tangential trajectory

%% get input arguments
% embedding vector
x = varargin{1};
N = size(x); % size of the embedding vector
if N(1) < N(2)
   error('Embedding dimension is larger than the length of the vector. Please check!')
end

% set threshold value
if nargin > 1
    if isa(varargin{2},'double')
        e = varargin{2};
    else
        warning('Threshold has to be numeric.')
    end
end

% set angle threshold value
if nargin > 2
    if isa(varargin{3},'double')
        w = varargin{3};
    else
        warning('Threshold has to be numeric.')
    end
end

% set time delay
if nargin > 3
    if isa(varargin{4},'double')
        tau = varargin{4};
    else
        warning('Time delay has to be numeric.')
    end
end


%% calculate distance matrix D and RP
% distance matrix using MATLAB's pdist function
d = pdist2(x,x);

% apply threshold and get the RP
r = d < e;


%% estimate the tangential vector
% two estimation variants:
% (1) use pre- and successor points
% (2) use reference point and another point in the future, tau time steps ahead
if 0
    % use pre- and successor points
    tangential_vec = zeros(size(x));
    tangential_vec(2:end-1,:) = x(3:end,:) - x(1:end-2,:);
else
    % use some delay (due to Horai et al, 2002)
    tangential_vec = zeros(size(x));
    tangential_vec(1:end-tau,:) = x(1+tau:end,:) - x(1:end-tau,:);
end

%% calculate dot product between the tangential vectors, if ~1 then parallel
sp = zeros(N(1),N(1));
for i = 1:N
   % dot product
   sp(i,:) = dot(tangential_vec, repmat(tangential_vec(i,:),N(1),1),2)';
   % normalize
   sp(i,:) = sp(i,:) ./ (vecnorm(tangential_vec,2,2) .* vecnorm(repmat(tangential_vec(i,:),N(1),1),2,2))';
end

% apply threshold to dot product matrix to check whether parallel (when dot product ~ 1)
iso_r = r .* (abs(sp-1) < w);

   
      

