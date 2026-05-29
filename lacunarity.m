function [L, Lnorm, Lshuffled, Lcomp] = lacunarity(X, boxSize)
%LACUNARITY Compute lacunarity measures for a binary or grayscale image.
%
%   [L, LNORM, LSHUFFLED] = LACUNARITY(X, BOXSIZE) computes the
%   lacunarity of the input matrix X for all box sizes specified in
%   BOXSIZE.
%
%   [L, LNORM, LSHUFFLED, LCOMP] additionally computes the lacunarity
%   of the complementary field (1 - X).
%
%   The function evaluates:
%
%     - Standard lacunarity
%     - Normalized lacunarity
%     - Lacunarity of a shuffled surrogate field
%     - Optional complementary lacunarity
%
%   Input:
%     X         : 2D numeric matrix (typically binary or grayscale)
%     boxSize   : Vector containing box sizes
%
%   Output:
%     L          : Standard lacunarity
%     Lnorm      : Normalized lacunarity
%     Lshuffled  : Lacunarity of shuffled surrogate data
%     Lcomp      : Lacunarity of complementary field (optional)
%
%   Example:
%       X = rand(2000) > 0.7;
%       X = checkerboard(200) > .5;
%       boxSize = 2:500;
%
%       [L, Lnorm, Lshuffled] = lacunarity(X, boxSize);
%
%       loglog(boxSize, L)
%       xlabel('Box size')
%       ylabel('Lacunarity')
%

% Copyright (c) 2026
% Potsdam Institute for Climate Impact Research, Germany
% Norbert Marwan
% https://www.pik-potsdam.de

L = zeros(length(boxSize), 1);
Lcomp = zeros(length(boxSize), 1);
Lshuffled = zeros(length(boxSize), 1);
Lnorm = zeros(length(boxSize), 1);

flag = 0;

if nargout > 3
    flag = 1;
end

% Create shuffled surrogate field
rng(42)

Xshuffle = reshape(X(randperm(numel(X))), size(X));

for i = 1:length(boxSize)

    % Box sums for original field
    n = boxsum(X, boxSize(i));

    m = mean(n(:));
    q = mean(n(:).^2);

    w2 = boxSize(i)^2;

    sigma2 = q - m^2;

    % Standard lacunarity
    L(i) = q / m^2;

    % Normalized lacunarity
    Lnorm(i) = 2 - ( ...
        m^2 / (m^2 + sigma2) + ...
        (w2 - m)^2 / ((w2 - m)^2 + sigma2) ...
    );

    % Complementary lacunarity
    if flag

        n = boxsum(1 - X, boxSize(i));

        Lcomp(i) = mean(n(:).^2) / mean(n(:))^2;

    end

    % Lacunarity of shuffled surrogate
    n = boxsum(Xshuffle, boxSize(i));

    Lshuffled(i) = mean(n(:).^2) / mean(n(:))^2;

end


function m = boxsum(X, boxSize)
%BOXSUM Compute local box sums of a 2D matrix.
%
%   M = BOXSUM(X, BOXSIZE) computes the summed values inside moving
%   square windows of size BOXSIZE x BOXSIZE over the input matrix X.
%
%   Depending on the box size, different algorithms are used:
%
%     - Small boxes : convolution-based method (conv2)
%     - Larger boxes: summed area table (integral image) method
%
%   For larger box sizes, the result may be subsampled internally to
%   improve performance.
%
%   Input:
%     X         : 2D numeric matrix
%     boxSize   : Size of the square moving window
%
%   Output:
%     m         : Matrix containing the local box sums
%
%   Example:
%       X = rand(100);
%       m = boxsum(X, 5);
%
%       imagesc(m)
%       colorbar
%       title('Local box sums')
%
%   See also CONV2, CUMSUM

% Copyright (c) 2026
% Potsdam Institute for Climate Impact Research, Germany
% Norbert Marwan
% https://www.pik-potsdam.de

method = 1;

if boxSize > 4
    step = round(boxSize / 1);
else
    step = 1;
end

if boxSize < 10
    method = 2;
end

switch(method)

    case 1
        %% Summed area table (integral image) approach

        % Compute integral image
        I = cumsum(cumsum(X, 1), 2);

        k = boxSize;
        s = step;

        % Valid starting indices
        rows = 1:s:(size(I,1) - boxSize + 1);
        cols = 1:s:(size(I,2) - boxSize + 1);

        % Prepare coordinates
        r1 = rows;
        r2 = rows + k - 1;

        c1 = cols;
        c2 = cols + k - 1;

        % Compute box sums
        m = I(r2, c2);

        % Subtract upper strip (for rows > 1)
        if any(r1 > 1)
            idx = r1 > 1;
            m(idx, :) = m(idx, :) - I(r1(idx)-1, c2);
        end

        % Subtract left strip (for cols > 1)
        if any(c1 > 1)
            idx = c1 > 1;
            m(:, idx) = m(:, idx) - I(r2, c1(idx)-1);
        end

        % Add overlapping corner region back
        if any(r1 > 1) && any(c1 > 1)
            ridx = r1 > 1;
            cidx = c1 > 1;
            m(ridx, cidx) = m(ridx, cidx) + ...
                I(r1(ridx)-1, c1(cidx)-1);
        end

    case 2
        %% Convolution-based approach

        kernel = ones(boxSize);

        sumX = conv2(X, kernel, 'valid');

        m = sumX(1:step:end, 1:step:end);

end
