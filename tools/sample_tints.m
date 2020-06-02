function t = sample_tints(tints, n)
%GETTRAINPOINTS Summary of this function goes here
%   Uniform sampling of n points across disjoint time intervals (tints)
    t_diff = diff(tints.epoch);
    t_diff = t_diff(1:2:end);
    t_diff_cumulative = cumsum(t_diff);
    t_diff_tot = double(t_diff_cumulative(end));
    [~, idx] = max(t_diff_cumulative' >= rand(n, 1)*t_diff_tot, [], 2);
    t_start = tints(2*idx - 1);
    t_diff = t_diff(idx);
    t = t_start + ((double(t_diff)*10^-9).*rand(n, 1));
end

