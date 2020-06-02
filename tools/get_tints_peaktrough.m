function tints = get_tints_peaktrough(data, data_mean, t)
%GET_TINTS_PEAKTROUGH Summary of this function goes here
%   Detailed explanation goes here
    dif = diff(data > data_mean);
    idx = find(dif)';
    idx = [1, reshape([idx; idx + 1], 1, []), length(data)];
    if data(1,:) > data_mean(1,:)
        idx = idx(3:end);
    end
    if data(end,:) > data_mean(end,:)
        idx = idx(1:end-2);
    end
    tints = t(idx);
end

