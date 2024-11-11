# Recurrence Plot & Quantification

[![Build Status](https://travis-ci.com/pucicu/rp.svg?branch=master)](https://travis-ci.com/pucicu/rp)
![commits](https://badgen.net/github/release/pucicu/rp)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/pucicu/rp/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/pucicu/rp)
![commits](https://badgen.net/github/license/pucicu/rp)



Simple MATLAB functions for calculating recurrence plots and recurrence quantification.

## Functions

### EMBED

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

```matlab
N = 300; % length of time series
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

```matlab
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

```matlab
N = 300; % length of time series
x = .9*sin((1:N)*2*pi/70); % exemplary time series
xVec = embed(x,2,17);
R = rp(xVec,.1);
Y = rqa(R);
```

--------------------------------------------------------------

### RP_ISO

Calculates the isodirectional recurrence plot

#### Syntax

`R=RP_ISO(X,E,W)` calculates the isodirectional recurrence plot `R`
from an embedding vector `X` and using the threshold `E` for the
vector distances and threshold `W` for the angle to be
considered as isodirectional.

`R=RP_ISO(X,E,W,TAU)` estimates tangential vector using time delay `TAU`.

#### Example

```matlab
[t x] = ode45('lorenz',[0 100],[-6.2 -10 14]);
[R1, SP, R0] = rp_iso(x(3000:5000,:),10,.2);

nexttile
imagesc(R0) % regular RP
axis square

nexttile
imagesc(R1) % isodirectional RP
axis square
```

--------------------------------------------------------------

### RP_PERP

Calculates the perpendicular recurrence plot

#### Syntax

`R=RP_PERP(X,E,W)` calculates the perpendicular recurrence plot `R`
from an embedding vector `X` and using the threshold `E` for the
vector distances and threshold `W` for the angle to be
considered as perpendicular.

`R=RP_PERP(X,E,W,TAU)` estimates tangential vector using time delay `TAU`
(works only if condition in line 95 is set to 0).

#### Example

```matlab
[t x] = ode45('lorenz',[0 100],[-6.2 -10 14]);
[R1, SP, R0] = rp_perp(x(3000:5000,:),10,.25);

nexttile
imagesc(R0) % regular RP
axis square

nexttile
imagesc(R1) % perpendicular RP
axis square
```

## Application

Part of this code was used in the study

* M. H. Trauth, A. Asrat, W. Duesing, V. Foerster, K. H. Kraemer, N. Marwan, M. A. Maslin, F. Schaebitz: _Classifying past climate change in the Chew Bahir basin, southern Ethiopia, using recurrence quantification analysis_, Climate Dynamics, 53(5), 2557–2572 (2019). DOI:[10.1007/s00382-019-04641-3](https://doi.org/10.1007/s00382-019-04641-3)

## How to cite

* Use the citation provided by Zenodo <https://zenodo.org/record/6325324>

## Background

read more about recurrence plot analysis at
<http://www.recurrence-plot.tk/>

## License

(see LICENSE file)

Copyright 2016-2020,
Potsdam Institute for Climate Impact Research (PIK),
Institute of Geosciences, University of Potsdam,
K. Hauke Kraemer, Norbert Marwan, Martin H. Trauth
