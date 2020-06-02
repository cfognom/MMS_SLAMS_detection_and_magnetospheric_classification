function tints = merge_tints(tints, lim)
%MERGE_TINTS Summary of this function goes here
%   Detailed explanation goes here
    if isempty(tints)
        tints = [];
        return
    else
        dur = diff(tints.epoch);
        between_dur = dur(2:2:end);
        keep = (between_dur > 1e9*lim)';
        keep = [1, reshape([keep; keep], 1, []), 1];
        tints = tints(keep == 1);
        % dur = dur(1:2:end);
        % dt = func(dur)';
        % dt = reshape([-dt; dt], 1, []);
        % tints = tints + dt';
        % tints = union_tints(tints);
    end
end

