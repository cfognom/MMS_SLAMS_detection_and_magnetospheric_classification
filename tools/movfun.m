function out = movfun(fun, data, dim, windowSize)
%MOVFUN Summary of this function goes here
%   Detailed explanation goes here
    if ~bitget(windowSize,1)
        error('Only odd window sizes are allowed.')
    end

    if dim == 2
        data = data';
    end

    b = (windowSize - 1)/2;
    s = size(data);
    n = s(1);
    hank = hankel(1 - b:n - b, n - b:n + b);
    out = zeros(n, 1);
    for i = 1:n
        idx = hank(i, :);
        idx = idx(idx >= 1 & idx <= n);
        widw = data(idx, :);
        out(i,:) = fun(widw);
    end

    if dim == 2
        out = out';
    end
end

