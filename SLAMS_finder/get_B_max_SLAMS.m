function [B_max, B_max_rel] = get_B_max_SLAMS(tints_SLAMS, b_abs, b_bg)
%GET_B_MAX_SLAMS Summary of this function goes here
%   Detailed explanation goes here
    if isempty(tints_SLAMS)
        B_max = [];
        B_max_rel = [];
        return
    end

    n_SLAMS = length(tints_SLAMS)/2;
    B_max = zeros(n_SLAMS, 1);
    B_max_rel = zeros(n_SLAMS, 1);
    for i = 1:n_SLAMS
        tint_SLAMS = select_tint(tints_SLAMS, i);
        b_data_SLAMS = b_abs.tlim(tint_SLAMS).data;
        b_rel_data_SLAMS = b_data_SLAMS./b_bg.tlim(tint_SLAMS).data;
        B_max(i) = max(b_data_SLAMS, [], 1);
        B_max_rel(i) = max(b_rel_data_SLAMS, [], 1);
    end
end

