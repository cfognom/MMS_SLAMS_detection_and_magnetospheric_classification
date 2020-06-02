function tints_uniq = labels2tints(ts_label, uniq)
%LABELS2TINTS Summary of this function goes here
%   Detailed explanation goes here
    y = ts_label.data';
    t = ts_label.time.epoch';
    n_t = length(t);
    n_uniq = length(uniq);
    id = find(diff(y));
    t_between = (t(id) + t(id + 1))/2;
    t_edges = [t(1), reshape([t_between; t_between], 1, []), t(n_t)];
    idx = [1, reshape([id; id + 1], 1, []), n_t];
    y_edges = y(idx);
    tints_uniq = cell(1, n_uniq);
    for i = 1:n_uniq
        tints_uniq{i} = EpochTT(int64(t_edges(y_edges == uniq(i))));
        % tints_uniq{i}
        % length(tints_uniq{i}) == length(unique(tints_uniq{i}.epoch))
    end
end

