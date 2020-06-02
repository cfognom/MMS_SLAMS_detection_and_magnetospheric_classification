function dur = get_duration_tints(tints)
%GET_DURATION_TINTS Summary of this function goes here
%   Detailed explanation goes here
    if isempty(tints)
        dur = 0;
        return
    end
    dif = diff(tints.epoch);
    dif = dif(1:2:end);
    dur = sum(dif);
end

