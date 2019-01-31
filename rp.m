function [y,P] = rp(varargin)
%
%   Version 1.4
%   new: - threshold-value in case of 'var'-setting now based on lower
%          quantile of the according distance distribution gained from the
%          adjacency-matrix
%
%
% Minimum input-arguments : 2
% Maximum input-arguments : 5
%
%    [R,DM] = RP(Y,E,thres-calc,norm,algorithm) 
%    
%    Calculates a recurrence plot R from an embedding vector Y and using 
%    the threshold 'E'. 'DM' is an optional output and is the adjacency- or
%    distancematrix.
%
% Input:
%
% 'Y'                       is a N-by-M matrix corresponding to N time points
%                           and M embedding dimensions.
%
% 'norm'                    norm for distance calculation in phasespace to
%                           'euc' (euclidic) or 'max' (maximum). Default is 
%                           max norm.
% 'algorithm'               specify the way of calculating the distance
%                           matrix here. You can choose from
%                           ['loops','vector','matlabvector']. Default is
%                           'vector'.
%
% 'threshold-calc'  specifies how the threshold epsilon will
% be calculated. There are three options. Set 'threshold-calc' to
%   - 'fix' The RP is computed under a fixed threshold epsilon specified by
%           input parameter 'E'.
%   - 'var' The RP is computed under a fixed threshold epsilon, which
%           corresponds to the lower 'E'-quantile (specified by input parameter
%           'E') of the distance distribution of all points in phasespace.
%   - 'fan' The RP is computed under a variable threshold epsilon using a
%           fixed amount of nearest neighbours in phasespace to compute the
%           epsilon-value for each point of the phasespace trajectory
%           individually.
% Default is 'fix'.  
%
%
%    Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         xVec = embed(x,2,17);
%         R = rp(xVec,.1);
%         imagesc(R)

% Copyright (c) 2016
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
% Modified by Hauke Krämer,Potsdam Institute for Climate Impact Research, 
% Germany http://www.pik-potsdam.de
% Institute of Earth and Environmental Science, University of Potsdam,
% Germany
%
% Contact: hkraemer@pik-potsdam.de, hkraemer@uni-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


%% check input
narginchk(1,5)
nargoutchk(0,2)

algoLib={'loops','vector','matlabvector'}; % the possible algorithms
try
    algorithm = varargin{5};
    if ~isa(algorithm,'char') || ~ismember(algorithm,algoLib)
        warning(['Specified algorithm should be one of the following possible values:',...
           10,sprintf('''%s'' ',algoLib{:})])
    end
catch
    algorithm = 'vector';
end

methLib={'euc','max'}; % the possible norms
try
    meth = varargin{4};
    if ~isa(meth,'char') || ~ismember(meth,methLib)
       warning(['Specified norm should be one of the following possible values:',...
           10,sprintf('''%s'' ',methLib{:})])
    end
catch
    meth = 'max';
end

thresLib={'fix','var','fan'}; % the possible ways of threshold computation
try
    thres = varargin{3};
    if ~isa(thres,'char') || ~ismember(thres,thresLib)
       warning(['Specified way of calculating threshold should be one of the following possible values:',...
                                10,sprintf('''%s'' ',thresLib{:})])
    end
catch
    thres = 'fix';
end

try
    e = varargin{2};
catch
    e = 1;
end

x = varargin{1};

N = size(x);
if N(1) < N(2)
   warning('Embedding dimension is larger than the length of the vector. Please check!')
end

%% init output

%% bind length of input vector in order to have constant iteration bounds while using parallel computing
M=N(1);
%% calculate distance matrix
switch algorithm 
    case 'loops'
         %% calculation with loops
         y = zeros(N(1),N(1));
         parfor i = 1:M
               for j = 1:M
                switch lower(meth)
                    case 'euc'
                         d = (x(i,:) - x(j,:)).^2;
                         y(i,j) = sqrt(sum(d));
                    otherwise
                         d = abs(x(i,:) - x(j,:));
                      y(i,j) = max(d);
                end
               end
         end
   case 'vector'
    
        %% calculation with vectorisation
        x1 = repmat(x,N(1),1);
        x2 = reshape(repmat(reshape(x,N(1)*N(2),1),1,N(1))',N(1)*N(1),N(2));
        switch lower(meth)
          case 'euc'
              d = (x1 - x2).^2;
              y = sqrt(sum(d,2));
          otherwise
              d = abs(x1 - x2);
              y = max(d,[],2);
        end
        y = reshape(y,N(1), N(1));   
        
    case 'matlabvector'
      %% calculation with matlab's vectorisation
      switch lower(meth)
          case 'euc'
              y = squareform(pdist(x));
          otherwise
              y = squareform(pdist(x,'chebychev'));
      end
end
P = y;
if strcmp(thres,'fix')
    % apply threshold
    y = double(y < e);
elseif strcmp(thres,'var')
    % apply threshold
    
    % lower (e*100)%-quantile of distance-distribution
    
    % transform adjacency-matrix into line-vector
    distance_vector = zeros(1,size(y,1)*size(y,2));
    start = 1;
    for i = 1:size(y,1)
        distance_vector(start:start + size(y,2)-1) = y(i,:);
        start = start + size(y,2);
    end
    % sort this line-vector
    distance_vector = sort(distance_vector);
    
    % get the lower 5%-quantile
    quantile_index = round(length(distance_vector)*e);
    epsilon = distance_vector(quantile_index);
    y = double(y < epsilon);
    

elseif strcmp(thres,'fan')
    % calculate threshold for each point 
    thresholds=zeros(M,M);
    for i = 1:M
        helper=sort(y(i,:));
          thresholds(i,:)=helper(round(e*M)+2);
    end
    % apply threshold(s)
    y=double(y<thresholds);
end
