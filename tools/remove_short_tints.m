function tints = remove_short_tints(tints, time_limit)
%REMOVE_SMALL_TINTS Summary of this function goes here
%   Detailed explanation goes here
    if isempty(tints)
        tints = [];
    else
        dif = diff(tints.epoch);
        dif = dif(1:2:end);
        idx = 2*find(dif > time_limit*10^9)';
        idx = reshape([idx - 1; idx], 1, []);
        tints = tints(idx);
    end
end

