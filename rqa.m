function y = rqa(varargin)
%RQA   Calculates recurrence quantification analysis
%    Q=RQA(R,L,T) calculates measures of recurrence 
%    quantification analysis for recurrence plot R using 
%    minimal line length L and a Theiler window T. 
%
%    Output:
%      Y(1) = RR     (recurrence rate)
%      Y(2) = DET    (determinism)
%      Y(3) = <L>    (mean diagonal line length)
%      Y(4) = Lmax   (maximal diagonal line length)
%      Y(5) = ENTR   (entropy of the diagonal line lengths)
%      Y(6) = LAM    (laminarity)
%      Y(7) = TT     (trapping time)
%      Y(8) = Vmax   (maximal vertical line length)
%      Y(9) = RTmax  (maximal white vertical line length)
%      Y(10) = T2     (recurrence time of 2nd type)
%      Y(11) = RTE    (recurrence time entropy, i.e., RPDE)
%      Y(12) = Clust  (clustering coefficient)
%      Y(13) = Trans  (transitivity)
% 
%    Reference:
%         Marwan, N., Romano, M. C., Thiel, M., Kurths, J. (2007).
%         Recurrence plots for the analysis of complex systems.
%         Physics Reports, 438, 237-329. 
%         Marwan, N., Donges, J. F., Zou, Y., Donner, R. V., 
%         Kurths, J. (2009). Complex network approach for recurrence
%         analysis of time series. Physics Letters A, 373, 4246-4254. 
%
%    Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         xVec = embed(x,2,17);
%         R = rp(xVec,.1);
%         Y = rqa(R);

% Copyright (c) 2016-2019
% Potsdam Institute for Climate Impact Research, Germany
% Institute of Geosciences, University of Potsdam, Germany
% Norbert Marwan, Hauke Kraemer
% http://www.pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%% check input
narginchk(1,3)
nargoutchk(0,1)

%% set default values for input parameters
theiler_window = 1; % Theiler window
l_min = 2; % minimal line length

%% get input arguments
% Theiler window
if nargin > 2
    theiler_window = varargin{3};
end

% minimal line length
if nargin > 1
    l_min = varargin{2};
end

% recurrence matrix
x = varargin{1};
N = size(x); % size of recurrence plot


%% apply Theiler window to recurrence plot
if theiler_window
   x_theiler = double(triu(x,theiler_window) + tril(x,-theiler_window));
else
   x_theiler = double(x);
end

% reduce the number of possible recurrence points by the Theiler window
N_all = N(1)*N(2); % all possible recurrence points
N_all = N_all - N(1) - 2*((theiler_window-1)*N(1) - sum(1:(theiler_window-1))); % reduced by Theiler window


%% calculation
y = zeros(13,1); % allocate result matrix

% recurrence rate
N_recpoints = sum(x_theiler(:)); % number of rec. points (in complete RP)
y(1) = N_recpoints/N_all; 


% histogram of diagonal lines (look at complete RP)
l_hist = zeros(1,N(1)); % allocate vector
for i = (1+theiler_window):N(1) % walk along the raws (upper triangle)
   cnt = 0;
   for j = 1:(N(2)-(i-1))
      if x_theiler(i+j-1,j) % are we on a rec. point? (walk along a diagonal)
         cnt = cnt+1; % count number of points on a diagonal line
      else % line has ended
         if cnt 
             l_hist(cnt) = l_hist(cnt) + 1; % store line length
         end
         cnt = 0; % set back to zero for a new line
      end
   end
   if cnt 
       l_hist(cnt) = l_hist(cnt) + 1;
   end
end

% 2nd triangle
for j = (1+theiler_window):N(2)  % walk along the columns (lower triangle)
   cnt = 0;
   for i = 1:(N(1)-(j-1))
      if x_theiler(i,j+i-1) % are we on a rec. point? (walk along a diagonal)
         cnt = cnt+1; % count number of points on a diagonal line
      else % line has ended
         if cnt 
             l_hist(cnt) = l_hist(cnt) + 1; % store line length
         end
         cnt = 0; % set back to zero for a new line
      end
   end
   if cnt
       l_hist(cnt) = l_hist(cnt) + 1; 
   end
end


% determinism
y(2) = sum(l_hist(l_min:N(1)) .* (l_min:N(1))) / N_recpoints;

% mean diagonal line length
if isnan(sum(l_hist(l_min:N(1)) .* (l_min:N(1))) / sum(l_hist(l_min:N(1))))
    y(3) = 0;
else
    y(3) = sum(l_hist(l_min:N(1)) .* (l_min:N(1))) / sum(l_hist(l_min:N(1)));
end

% maximal line length
if any(l_hist)
   y(4) = find(l_hist,1,'last');
end

% line length entropy
l_classes = sum(l_hist~=0); % number of occupied bins (for normalization of entropy)
l_prob = l_hist/sum(l_hist); % get probability distribution from histogram
ent_Sum = (l_prob .* log(l_prob));
if l_classes > 1
    y(5) = -nansum(ent_Sum)/log(l_classes);
else
    y(5) = -nansum(ent_Sum);
end


% histogram of vertical lines
v_hist = zeros(1,N(1)); % allocate vector
for i = 1:N(1) % walk along the columns
   cnt = 0;
   for j = 1:N(2)
      if x_theiler(j,i) % are we on a rec. point? (walk along a column)
         cnt = cnt+1; % count number of points on a vertical line
      else % line has ended
         if cnt
             v_hist(cnt) = v_hist(cnt) + 1; % store line length
         end
         cnt = 0; % set back to zero for a new line
      end
   end
   if cnt
       v_hist(cnt) = v_hist(cnt) + 1;
   end
end

% laminarity
y(6) = sum(v_hist(l_min:N(1)) .* (l_min:N(1))) / N_recpoints;

% mean vertical line length (trapping time)
if isnan(sum(v_hist(l_min:N(1)) .* (l_min:N(1))) / sum(v_hist(l_min:N(1))))
    y(7) = 0;
else
    y(7) = sum(v_hist(l_min:N(1)) .* (l_min:N(1))) / sum(v_hist(l_min:N(1)));
end

% maximal vertical length
if any(v_hist)
   y(8) = find(v_hist,1,'last');
end


% recurrence times ("white" vertical lines)
rt_hist = zeros(1,N(1)); % allocate vector
for i = 1:N(1)
   cnt = 0;
   
   % boolean variable to avoid counting white lines at the edges of RP
   first_flag = false;
   
   for j = 1:N(2)
      if ~x(j,i) % are we on a white line?
         if first_flag % line does not cross the RP's edges
            cnt = cnt + 1; % count number of points along the vertical line
         end
      else % we meet a recurrence point
         first_flag = true; % we are for sure within the RP
         if cnt
             rt_hist(cnt) = rt_hist(cnt) + 1; % store line length
         end
         cnt = 0;
      end
   end
end

% maximal white vertical line length
if any(rt_hist)
    y(9) = find(rt_hist,1,'last');
end

% recurrence time
y(10) = sum(rt_hist .* (1:N(1))) / sum(rt_hist);
if isnan(y(10))
    y(10)=0;
end

% recurrence time entropy
rt_classes = sum(rt_hist~=0); % number of occupied bins (for normalization of entropy)
rt_prob = rt_hist/sum(rt_hist); % get probability distribution from histogram
ent_Sum = (rt_prob .* log(rt_prob));
if rt_classes > 1
    y(11) = -nansum(ent_Sum)/log(rt_classes);
else
    y(11) = -nansum(ent_Sum);
end

% clustering
kv = sum(x_theiler,1); % degree of nodes
y(12) = nanmean(diag(x_theiler*x_theiler*x_theiler)' ./ (kv .* (kv-1)));

% transitivity
denom = sum(sum(x_theiler * x_theiler));
y(13) = trace(x_theiler*x_theiler*x_theiler)/denom;
