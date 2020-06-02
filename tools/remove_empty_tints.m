function tints = remove_empty_tints(tints)
%REMOVE_EMPTY_TINTS Summary of this function goes here
%   Detailed explanation goes here
    if isempty(tints)
        tints = [];
        return
    end
    dif = diff(tints.epoch);
    dif = dif(1:2:end);
    idx = 2*find(dif);
    idx = sort([idx - 1; idx]);
    tints = tints(idx);
    if isempty(tints)
        tints = [];
    end
end

