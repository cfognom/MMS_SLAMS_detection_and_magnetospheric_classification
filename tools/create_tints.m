function tints = create_tints(starts, stops)
%CREATE_TINTS Summary of this function goes here
%   Detailed explanation goes here
    n_tints = length(starts);

    idx = zeros(1, 2*n_tints);
    idx(1:2:end) = 1:n_tints;
    idx(2:2:end) = n_tints + 1:2*n_tints;

    tints = [starts; stops];
    tints = tints(idx);
end

