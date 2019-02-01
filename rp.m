function [r,d] = rp(varargin)
%RP   Calculate a recurrence plot
%    R=RP(X,E,THRESH,NORM,ALG) calculates the recurrence plot R 
%    from an embedding vector X and using the threshold E.
%    X is a N-by-M matrix corresponding to N time points
%    and M embedding dimensions.
%
%    [R,D]=RP(...) outputs the recurrence plot R and the 
%    underlying distance matrix D.
%
%    Optional arguments:
%          NORM - is a string setting the norm for distance 
%                calculation in phasespace. Can be 'euc' 
%                for euclidian norm (default) or 'max' 
%                for maximum norm.
%    ALGORITHM - is a string specifying the algorithm of 
%                calculating the distance matrix. Can be
%                'loops', 'vector' (default), or 
%                'matlabvector'.
%       THRESH - is a string that specifies how the threshold 
%                epsilon will be calculated. With 'fix' (default)
%                the RP is computed with a fixed threshold 
%                epsilon specified by the input parameter E.
%                With 'var' the RP is computed with a fixed 
%                threshold epsilon, which corresponds to the 
%                lower 'E'-quantile (specified by E) of the 
%                distance distribution of all points in 
%                phasespace. With 'fan' the RP is computed with
%                a variable threshold resulting in a fixed amount 
%                of nearest neighbours in phasespace.
%
%    Reference:
%         Marwan, N., Romano, M. C., Thiel, M., Kurths, J. (2007).
%         Recurrence plots for the analysis of complex systems.
%         Physics Reports, 438, 237-329. 
%         Kraemer, K. H., Donner, R. V., Heitzig, J., & Marwan, N. 
%         (2018). Recurrence threshold selection for obtaining robust
%         recurrence characteristics in different embedding dimensions.
%         Chaos, 28, 085720.
%
%    Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         xVec = embed(x,2,17); % embed into 2 dimensions using delay 17
%         R = rp(xVec,.1,'fix','max'); % calculate RP using maximum norm and fixed threshold
%         imagesc(R)

% Copyright (c) 2016-2019
% Potsdam Institute for Climate Impact Research, Germany
% Institute of Geo Sciences, University of Potsdam, Germany
% Norbert Marwan, Hauke Krämer
% http://www.pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


%% check input
narginchk(1,5)
nargoutchk(0,2)

%% set default values for input parameters
algoLib={'loops','vector','matlabvector'}; % the possible algorithms
methLib={'euc','max'}; % the possible norms
thresLib={'fix','var','fan'}; % the possible ways of threshold computation

algorithm = 'vector';
meth = 'max'; % norm 
thres = 'fix'; % threshold algorithm
e = 1; % recurrence threshold


%% get input arguments
% set the algorithm of computing the RP
if nargin > 4
    if isa(varargin{5},'char') & ismember(varargin{5},algoLib)
        algorithm = varargin{5};
    else
        warning(['Specified algorithm should be one of the following possible values:',...
           10,sprintf('''%s'' ',algoLib{:})])
    end
end

% set norm
if nargin > 3
    if isa(varargin{4},'char') & ismember(varargin{4},methLib)
        meth = lower(varargin{4});
    else
        warning(['Specified norm should be one of the following possible values:',...
           10,sprintf('''%s'' ',methLib{:})])
    end
end

% set method for threshold calculation
if nargin > 2
    if isa(varargin{3},'char') & ismember(varargin{3},thresLib)
        thres = lower(varargin{3});
    else
        warning(['Specified way of calculating threshold should be one of the following possible values:',...
            10,sprintf('''%s'' ',thresLib{:})])
    end
end

% set threshold value
if nargin > 1
    if isa(varargin{2},'double')
        e = varargin{2};
    else
        warning('Threshold has to be numeric.')
    end
end

% embedding vector
x = varargin{1};
N = size(x); % size of the embedding vector
if N(1) < N(2)
   error('Embedding dimension is larger than the length of the vector. Please check!')
end


%% calculate distance matrix D

switch algorithm 
    % calculation using loops
    case 'loops'
         d = zeros(N(1),N(1)); % allocate matrix
         for i = 1:N(1)
             for j = 1:N(1)
                 switch meth % select norm
                     case 'euc'
                         % euclidean distance between two phase space vectors
                         s = (x(i,:) - x(j,:)).^2;  
                         d(i,j) = sqrt(sum(s));
                     otherwise
                         % max-norm distance between two phase space vectors
                         s = abs(x(i,:) - x(j,:));
                         d(i,j) = max(s);
                 end
             end
         end
   % calculation using full vectorisation
   case 'vector'
        % repeat vector entries to have all pair-wise combinations
        x1 = repmat(x,N(1),1);
        x2 = reshape(repmat(reshape(x,N(1)*N(2),1),1,N(1))',N(1)*N(1),N(2));
        switch meth
            case 'euc'
                % euclidean distance between two phase space vectors
                s = (x1 - x2).^2;
                d = sqrt(sum(s,2));
            otherwise
                % max-norm distance between two phase space vectors
                s = abs(x1 - x2);
                d = max(s,[],2);
        end
        d = reshape(d,N(1), N(1)); % reshape to have a square matrix
        
    % calculation using MATLAB's pdist function
    case 'matlabvector'
      switch meth
          case 'euc'
              % euclidean distance between two phase space vectors
              d = squareform(pdist(x));
          otherwise
              % max-norm distance between two phase space vectors
              d = squareform(pdist(x,'chebychev'));
      end
end


%% apply threshold to get the final RP


% calculate the threshold
if strcmp(thres,'fix')
    % simply apply threshold to the distance matrix
    r = double(d < e);

elseif strcmp(thres,'var')
    % threshold based on lower (e*100)%-quantile of distance-distribution
        
    % get the lower e%-quantile
    epsilon = quantile(d(:),e);
    r = double(d < epsilon);

elseif strcmp(thres,'fan')
    % variable threshold for each point in order to get fixed
    % number of nearest neighbours
    q = quantile(d,e); % distance that corresponds to the fraction e of rec. points per column
    thresholds = repmat(q,N(1),1); % q has to be applied for each row in d
    % apply individual thresholds
    r = double(d<thresholds);
end

