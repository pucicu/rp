function y = rqa(varargin)
%
%   version 1.4
%   new:        - normalized Entropy-measurements corrected
%
%
%    Calculate recurrence quantification analysis
%    Y=RQA(X,L,T) calculates measures of recurrence 
%    quantification analysis using minimal line length L
%    and a Theiler window T. The reus
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
%
%    Example:
%         N = 300; % length of time series
%         x = .9*sin((1:N)*2*pi/70); % exemplary time series
%         xVec = embed(x,2,17);
%         R = rp(xVec,.1);
%         Y = rqa(R);

% Copyright (c) 2016
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Modified by Hauke Krämer, Potsdam Institute for Climate Impact Research,
% Germany
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%% check input
narginchk(1,3)
nargoutchk(0,1)

try
    theiler_window = varargin{3};
catch
    theiler_window = 1;
end

try
    l_min = varargin{2};
catch
    l_min = 2;
end
x = varargin{1};

% size of recurrence plot
N = size(x);

%% calculation
y = zeros(13,1);

% applt Theiler window to recurrence plot
if theiler_window
   x_theiler = double(triu(x,theiler_window) + tril(x,-theiler_window));
else
   x_theiler = double(x);
end



% reduce the number of possible states by the Theiler window
N_all = N(1)*N(2);
N_all = N_all - N(1) - 2*((theiler_window-1)*N(1) - sum(1:(theiler_window-1)));

% recurrence rate
N_recpoints = sum(x_theiler(:)); % number of rec. points (in complete RP)
y(1) = N_recpoints/N_all; 

% histogram of diagobal lines (look at entire RP)
l_hist = zeros(1,N(1));
for i = (1+theiler_window):N(1)
   cnt = 0;
   for j = 1:(N(2)-(i-1))
      if x_theiler(i+j-1,j)
         cnt = cnt+1; % count number of points on a diagonal line
      else
         if cnt 
             l_hist(cnt) = l_hist(cnt) + 1; 
         end
         cnt = 0;
      end
   end
   if cnt 
       l_hist(cnt) = l_hist(cnt) + 1;
   end
end

for j = (1+theiler_window):N(2)
   cnt = 0;
   for i = 1:(N(1)-(j-1))
      if x_theiler(i,j+i-1)
         cnt = cnt+1; % count number of points on a diagonal line
      else
         if cnt 
             l_hist(cnt) = l_hist(cnt) + 1;
         end
         cnt = 0;
      end
   end
   if cnt
       l_hist(cnt) = l_hist(cnt) + 1; 
   end
end
% check how many classes occupied in diagonal line histogram
l_classes = 0;
for i = 1:length(l_hist)
    if l_hist(i)~=0
        l_classes = l_classes + 1;
    end
end

% Determinism
N_recpoints2 = sum(l_hist(1:N(1)) .* (1:N(1)));  % number of rec. points (in one triangle of RP)
y(2) = sum(l_hist(l_min:N(1)) .* (l_min:N(1))) / N_recpoints2;

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
l_prob = l_hist/sum(l_hist); % get probability distribution from histogram
ent_Sum = (l_prob .* log(l_prob));
y(5) = -sum(ent_Sum(~isnan(ent_Sum)))/log(l_classes);

% laminarity
% histogram of vertical lines
v_hist = zeros(1,N(1));
for i = 1:N(1)
   cnt = 0;
   for j = 1:N(2)
      if x_theiler(i,j)
         cnt = cnt+1; % count number of points on a vertical line
      else
         if cnt
             v_hist(cnt) = v_hist(cnt) + 1; 
         end
         cnt = 0;
      end
   end
   if cnt
       v_hist(cnt) = v_hist(cnt) + 1;
   end
end

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

rt_hist = zeros(1,N(1));
for i = 1:N(1)
   cnt = 0;
   
   % boolean variable to avoid counting white lines at the edges of RP
   first_flag = false;
   
   for j = 1:N(2)
      if ~x(i,j)
         if first_flag == true
            cnt = cnt + 1;
         end
      else
         first_flag = true;
         if cnt
             rt_hist(cnt) = rt_hist(cnt) + 1; 
         end
         cnt = 0;
      end
   end
end

% check how many classes occupied in white line histogram
rt_classes = 0;
for i = 1:length(rt_hist)
    if rt_hist(i)~=0
        rt_classes = rt_classes + 1;
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
rt_prob = rt_hist/sum(rt_hist); % get probability distribution from histogram
ent_Sum = (rt_prob .* log(rt_prob));
y(11) = -sum(ent_Sum(~isnan(ent_Sum)))/log(rt_classes);

% clustering
kv = sum(x_theiler,1); % degree of nodes
y(12) = nanmean(diag(x_theiler*x_theiler*x_theiler)' ./ (kv .* (kv-1)));

% transitivity
denom = sum(sum(x_theiler * x_theiler));
y(13) = trace(x_theiler*x_theiler*x_theiler)/denom;
