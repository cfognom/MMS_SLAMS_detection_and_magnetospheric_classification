function disc = discretize_tints(tints_bins, tints)
%DISCRETIZE_TINTS Summary of this function goes here
%   Detailed explanation goes here
    n_bins = length(tints_bins)/2;
    disc = zeros(1, n_bins);
    for i = 1:n_bins
        tint_bin = select_tint(tints_bins, i);
        tints_inside_bin = intersect_tints(tint_bin, tints);
        bin_dur = get_duration_tints(tint_bin);
        inside_bin_dur = get_duration_tints(tints_inside_bin);
        inside_percentage = inside_bin_dur/bin_dur;
        if inside_percentage > 0.5
            disc(i) = 1;
        end
    end
end

