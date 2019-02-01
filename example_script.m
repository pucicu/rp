%% MATLAB Script to RPs/RQA for Synthetic Data
% MATLAB script to reproduce the recurrence plots (RPs) and the recurrence
% quantification analysis (RQA) measures for synthetic data representing
% common types of dynamic behavior. This script creates Figure 3 of the
% paper:
%
% Trauth, M.H., Asrat, A., Duesing, W., Foerster, V., Kraemer, K.H.,
% Marwan, N., Maslin, M.A., Schaebitz, F. (2019) Classifying past climate
% change in the Chew Bahir basin, southern Ethiopia, using recurrence
% quantification analysis. Climate Dynamics, Springer Verlag GmbH Germany,
% https://doi.org/10.1007/s00382-019-04641-3
%
% and can be used to perform the RP/RQA of the Chew Bahir data presented in
% the paper.

% Copyright (c) 2016-2019
% Potsdam Institute for Climate Impact Research, Germany
% Institute of Geosciences, University of Potsdam, Germany
% Norbert Marwan, Hauke Kraemer, Martin H. Trauth
% http://www.pik-potsdam.de and http://www.geo.uni-potsdam.de
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or any later
% version.

%%
% Clear workspace, clear Command Window, close all Figure windows.
clear, clc, close all

%%
% Select RQA measures and create RQA labels.
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
RQA_Legend = [
        'RR   '
        'DET  '
        '<L>  '
        'LMAX '
        'ENTR '
        'LAM  '
        'TT   '
        'VMAX '
        'RTMAX'
        'T2   '
        'RTE  '
        'CLUST'
        'TRANS'];
RQA_Select = [1 2 13];
RQA_Labels = ['A','B','C','D','E','F'];

%%
% Define the method of threshold calculation for RQA. Calculate recurrence
% plot using the recurrence threshold e and choose threshold calculation
% parameter. Set parameter to define the threshold calculation for the
% recurrence plot estimation in the next step. There are three options to
% choose from:
%
%   - 'fix' The RP is computed under a fixed threshold epsilon specified by
%           input parameter "e".
%   - 'var' The RP is computed under a fixed threshold epsilon, which
%           corresponds to the lower 'E'-quantile (specified by input
%           parameter 'e') of the distance distribution of all points in
%           the phase space.
%   - 'fan' The RP is computed under a variable threshold epsilon using a
%           fixed amount of nearest neighbours in phasespace to compute the
%           epsilon-value for each point of the phasespace trajectory
%           individually.
threshold_calculation = 'var';

%%
% Choose norm from ['euc','max'], define Theiler window and minimal line
% length.
norm = 'euc';
theiler = 1;
l_min = 4;

%%
% Create synthetic data, embedding parameters, threshold and window.
t = -2000 : 2 : 0;

% Example 1: Normally-distributed (Gaussian) noise.
m = 1;
rng(0), y = randn(size(t));  
y = y/std(y);  % normalize time series
xx(m,:) = y;
mm(m) = 3;
tautau(m) = 1;
if strcmp(threshold_calculation,'fan')||strcmp(threshold_calculation,'var')
    ee(m) = 0.07;
else
    ee(m) = 1;
end   
ww(m) = 100;
wsws(m) = 10;

% Example 2: Composite signal comprising two sine waves and a positive
% trend in the mean.
m = 2;
y = sin(2*pi*t/200) + sin(2*pi*t/50);

counter = 1;
for i = 1:length(t)
    if i >=300 && i <= 700
        y(i) = y(i) + counter*0.0003*(i-349);
        base = counter*0.0003*(i-349);
        if counter < 50
            counter = counter + 1;
        end
    elseif i > 700
        y(i) = y(i) + base;
    end

end
y = y/std(y); % normalize time series
xx(m,:) = y;
mm(m) = 4;
tautau(m) = 10;
if strcmp(threshold_calculation,'fan')||strcmp(threshold_calculation,'var')
    ee(m) = 0.07;
else
    ee(m) = 1;
end 
ww(m) = 100;
wsws(m) = 10;

% Example 3: Composite signal comprising a sine wave and Gaussian noise
% with decreasing signal-to-noise ratio from left to right.
m = 3;
y = sin(2*pi*t/200) + fliplr(t).*randn(size(y))/2000;

y = y/std(y);
xx(m,:) = y;
mm(m)     = 4;
tautau(m) = 25;
if strcmp(threshold_calculation,'fan')||strcmp(threshold_calculation,'var')
    ee(m) = 0.07;
else
    ee(m) = 1;
end 
ww(m) = 100;
wsws(m) = 10;

% Example 4: Composite signal comprising two sine waves and a trend in the
% frequencies.
m = 4;
cnt=1;
for i = 0:2:2000
    y(cnt) = sin(2*pi*i^(1.35)/1500) + sin(2*pi*i^(1.35)/400);
    cnt=cnt+1;
end
y = y/std(y);
xx(m,:) = y';
mm(m)     = 5;
tautau(m) = 5;
if strcmp(threshold_calculation,'fan')||strcmp(threshold_calculation,'var')
    ee(m) = 0.07;
else
    ee(m) = 1;
end 
ww(m) = 100;
wsws(m) = 10;

% Example 5: Abrupt transition in the amplitudes of two sine waves.
m = 5;
counter = 1;
counter2= 1;
counter3= 0;
counter4= 0;
for i = 1:length(t)

    if i > 350 && i < 441
        y(i) = sin(2*pi*t(i)/300) + (1-(counter*0.0111))*sin(2*pi*t(i)/50);
        counter = counter + 1;
    elseif i > 440 && i < 501
        y(i) = (1-counter2*0.018)*sin(2*pi*t(i)/(300-(counter3*4)));
        if counter2<51
            counter2 = counter2+1;
        end
        
        counter3 = counter3 +1;
    elseif i > 500 && i < 551
        y(i) = (0.1+counter4*0.018)*sin(2*pi*t(i)/60);
        counter4 = counter4 +1;
    else
        y(i) = sin(2*pi*t(i)/300) + sin(2*pi*t(i)/50);
    end
end               
y(551:1001) = sin(2*pi*t(551:end)/60);
y(1:500) = y(1:500)/std(y(1:500));
y(501:end) = y(501:end)/std(y(501:end));
xx(m,:) = y;
mm(m) = 3;
tautau(m) = 7;
if strcmp(threshold_calculation,'fan')||strcmp(threshold_calculation,'var')
    ee(m) = 0.07;
else
    ee(m) = 1;
end 
ww(m) = 100;
wsws(m) = 10;

% Example 6: Normally-distributed (Gaussian) noise with a stepwise
% transition in the mean and a change in the autocorrelation prior to this
% transition.
m = 6;
rng(20);
y = [zeros(500,1);ones(501,1)];
y = y + .3*randn(1001,1);
for i = 400:500
   y(i) = .5 * y(i-1) + .5 * y(i);
end
y = y/std(y);
xx(m,:) = y';
mm(m) = 5;
tautau(m) = 2;
if strcmp(threshold_calculation,'fan')||strcmp(threshold_calculation,'var')
    ee(m) = 0.07;
else
    ee(m) = 1;
end 
ww(m) = 100;
wsws(m) = 10;

%%
% Loop examples x.
close all
for k = 1 : size(xx,1)
    clear x m tau e w ws
    x = xx(k,:);

    % Create embedding vector, define embedding dimension m and delay tau. The 
    % embedding delay must be even for correct phase corrections.
    m = mm(k);
    tau = tautau(k);

    timespan_diff = tau*(m-1);
    xVec = embed(x,m,tau);

    % Use initial value of threshold e.
    e = ee(k);

    % Calculate RP.
    R = rp(xVec,e,threshold_calculation,norm);

    % Calculate RQA measures on complete time series.
    r = rqa(R,l_min,theiler);

    % Calculate RQA measures in moving windows.
    clear r_win r_win_w r_win_e

    w = ww(k);     % Window size.
    ws = wsws(k);  % Windowing moving by lag ws, ws>1.
    r_win=zeros(13,ceil((length(R)-w)/ws));  % Preallocate memory.
    cnt = 1;       % Counter.
    for i = 1:ws:(length(R)-w)
       r_win(:,cnt) = rqa(R(i:(i+w),i:(i+w)),l_min,theiler);
       cnt = cnt+1;
    end

    % Phase correction for windowed measures by shifting the measures by half
    % the window size and sampled at the resolution ws. To draw a continuous
    % line the phase corrected measures are interpolated linearly using
    % fillmissing.
    r_win_w(1:size(r_win,1),1:size(R,1)) = NaN;
    r_win_w(:,1+w/2:ws:size(R,1)-w/2) = r_win;
    for i = 1 : 13
        r_win_w(i,1+w/2:size(R,1)-w/2) = ...
        fillmissing(r_win_w(i,1+w/2:size(R,1)-w/2),'linear');
    end

    % Phase correction for embedding delay by shifting the measures by half the
    % embedding delay. To draw a continuous line the phase corrected measures
    % are interpolated linearly using fillmissing.
    r_win_e(1:size(r_win,1),1:size(x,2)) = NaN;
    r_win_e(:,1+round(timespan_diff/2):size(x,2)- ...
        floor(timespan_diff/2)) = r_win_w;

    RR = NaN(size(x,2),size(x,2));
    RR(1+round(timespan_diff/2):size(x,2)-floor(timespan_diff/2),...
        1+round(timespan_diff/2):size(x,2)-floor(timespan_diff/2)) = R;

    % Display time series, recurrence plot and RQA measures.

    % Figure.
    figure('Position',[(k-1)*400 600 400 600],'Color',[1 1 1])

    % Time series.
    h(1) = axes('Position',[0.1 0.79 0.8 0.15],...
         'XLim',[min(t) max(t)],...
         'LineWidth',0.75,...
         'Box','On',...
         'XTickLabel',''); hold on
    line(t,x,...
         'LineWidth',1,...
         'Color',[0 0 0])

    % Recurrence plot.
    h(2) = axes('Position',[0.1 0.11 0.8 0.8],...
         'XLim',[min(t) max(t)],...
         'YLim',[min(t) max(t)],...
         'LineWidth',0.75,...
         'Box','On',...
         'XTick',[],...
         'YTick',[]); hold on
    axis square xy
    imagesc(t,t,RR)
    axes('Position',[0.1 0.11 0.8 0.8],...
         'XLim',[min(t) max(t)],...
         'YLim',[min(t) max(t)],...
         'LineWidth',0.75,...
         'Box','On',...
         'Color','none',...
         'XTick',[],...
         'YTick',[]), hold on
    axis square xy
    colormap([1 1 1; 0 0 0])

    % Parameters.
    str1 = ['m    = ',num2str(m)];
    str2 = ['tau  = ',num2str(tau)];
    str3 = ['e    = ',num2str(e),' (adaptive)'];
    str4 = ['w    = ',num2str(w)];
    str5 = ['ws   = ',num2str(ws)];
    str6 = ['norm = ',norm];
    str7 = ['thei = ',num2str(theiler)];
    str8 = ['lmin = ',num2str(l_min)];
    str = {str1,str2,str3,str4,str5,str6,str7,str8};
    text(min(t)+(max(t)-min(t))/20,...
         max(t)-(max(t)-min(t))/20,str,...
        'BackgroundColor',[1 1 1],...
        'EdgeColor',[1 1 1],...
        'VerticalAlignment','top')

    % RQA - Recurrence rate (blue) and determinism (red).
    h(3) = axes('Position',[0.1 0.08 0.8 0.15],...
        'XLim',[min(t) max(t)],...
        'YLim',[0 1],...
        'LineWidth',0.75,...
        'XGrid','On',...
        'Box','On'); hold on
    line(h(3),t,r_win_e(RQA_Select(1),:),...
        'LineWidth',1,...
        'Color',[0 0.4453 0.7383])
    line(h(3),t,r_win_e(RQA_Select(2),:),...
        'LineWidth',1,...
        'Color',[0.8477 0.3242 0.0977])
    legend({RQA_Legend(RQA_Select(1),:),...
           RQA_Legend(RQA_Select(2),:)},...
           'Box','Off')
    xlabel('Time')

    % Legend.
    legend({RQA_Legend(RQA_Select(1),:),...
           RQA_Legend(RQA_Select(2),:)},...
           'Box','Off')

    linkaxes(h,'x')

    % Print.
    printname = ['trauth_figure_3',RQA_Labels(k),'.png'];
    print(printname,'-dpng','-r300')

end


