# Recurrence plot & Quantification #

Simple MATLAB functions for calculating recurrence plots and recurrence quantification.


## Functions

### embed

Creates embedding vector using time delay embedding.

#### Syntax

`Y=EMBED(X,M,T)` creates the embedding vector `Y` from the time
series `X` using a time delay embedding with dimension `M` and
delay `T`. The resulting embedding vector has length `N-T*(M-1)`,
where `N` is the length of the original time series.

#### Reference

* Packard, N. H., Crutchfield, J. P., Farmer, J. D.,
  Shaw, R. S. (1980). Geometry from a time series.
  Physical Review Letters 45, 712-716.
          
#### Example         
```N = 300; % length of time series
x = .9*sin((1:N)*2*pi/70); % exemplary time series
y = embed(x,2,17); % embed into 2 dimensions using delay 17
plot(y(:,1),y(:,2))
```

--------------------------------------------------------------
### RP

Calculates a recurrence plot.


#### Syntax

`R=RP(X,E,THRESH,NORM,ALG)` calculates the recurrence plot `R`
from an embedding vector `X` and using the threshold `E`.
`X` is a `N`-by-`M` matrix corresponding to `N` time points
and `M` embedding dimensions.

`[R,D]=RP(...)` outputs the recurrence plot `R` and the
underlying distance matrix `D`.

**Optional arguments:**

`NORM` - is a string setting the norm for distance
                 calculation in phasespace. Can be `'euc'`
                 for euclidian norm (default) or `'max'`
                 for maximum norm.
                 
`ALG` - is a string specifying the algorithm of
                 calculating the distance matrix. Can be
                 `'loops'`, `'vector'` (default), or
                 `'matlabvector'`.

`THRESH` - is a string that specifies how the threshold
                 epsilon will be calculated. With `'fix'` (default)
                 the RP is computed with a fixed threshold
                 epsilon specified by the input parameter `E`.
                 With `'var'` the RP is computed with a fixed
                 threshold epsilon, which corresponds to the
                 lower '`E`'-quantile (specified by `E`) of the
                 distance distribution of all points in
                 phasespace. With `'fan'` the RP is computed with
                 a variable threshold resulting in a fixed amount
                 of nearest neighbours in phasespace, specified
%                by the fraction `E` of recurrence points

#### Reference
* Marwan, N., Romano, M. C., Thiel, M., Kurths, J. (2007).
  Recurrence plots for the analysis of complex systems.
  Physics Reports, 438, 237-329.
* Kraemer, K. H., Donner, R. V., Heitzig, J., & Marwan, N.
  (2018). Recurrence threshold selection for obtaining robust
  recurrence characteristics in different embedding dimensions.
  Chaos, 28, 085720.


#### Example
```
N = 300; % length of time series
x = .9*sin((1:N)*2*pi/70); % exemplary time series
xVec = embed(x,2,17); % embed into 2 dimensions using delay 17
R = rp(xVec,.1,'fix','max'); % calculate RP using maximum norm and fixed threshold
imagesc(R)
```      

--------------------------------------------------------------

### RQA

Calculates recurrence quantification analysis.

#### Syntax

`Q=RQA(R,L,T)` calculates measures of recurrence
quantification analysis for recurrence plot `R` using
minimal line length `L` and a Theiler window `T.

**Output:**

* `Y(1) = RR`     (recurrence rate)
* `Y(2) = DET`    (determinism)
* `Y(3) = <L>`    (mean diagonal line length)
* `Y(4) = Lmax`   (maximal diagonal line length)
* `Y(5) = ENTR`   (entropy of the diagonal line lengths)
* `Y(6) = LAM`    (laminarity)
* `Y(7) = TT`     (trapping time)
* `Y(8) = Vmax`   (maximal vertical line length)
* `Y(9) = RTmax` (maximal white vertical line length)
* `Y(10) = T2`     (recurrence time of 2nd type)
* `Y(11) = RTE`    (recurrence time entropy, i.e., RPDE)
* `Y(12) = Clust`  (clustering coefficient)
* `Y(13) = Trans`  (transitivity)

#### Reference
* Marwan, N., Romano, M. C., Thiel, M., Kurths, J. (2007).
  Recurrence plots for the analysis of complex systems.
  Physics Reports, 438, 237-329.
* Marwan, N., Donges, J. F., Zou, Y., Donner, R. V.,
  Kurths, J. (2009). Complex network approach for recurrence
  analysis of time series. Physics Letters A, 373, 4246-4254.

#### Example
```
N = 300; % length of time series
x = .9*sin((1:N)*2*pi/70); % exemplary time series
xVec = embed(x,2,17);
R = rp(xVec,.1);
Y = rqa(R);
```     

## Application

Trauth et al, Classifying past climate change in the Chew Bahir basin, southern Ethiopia, using recurrence quantification analysis, *Climate Dynamics*, in press, 2019
          
## How to cite

* Zenodo?

## Background

http://www.recurrence-plot.tk/