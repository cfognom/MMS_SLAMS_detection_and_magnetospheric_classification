function [n_tints, starts, stops] = unpack_tints(tints)
%TINT_SPLIT Summary of this function goes here
%   Detailed explanation goes here
    n_t = length(tints);
    starts = tints(1:2:n_t);
    stops = tints(2:2:n_t);
    n_tints = n_t/2;
end

