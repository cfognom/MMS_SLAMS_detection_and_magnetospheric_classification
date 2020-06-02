function tint = select_tint(tints, index)
%SELECT_TINT Summary of this function goes here
%   Detailed explanation goes here
    index = 2*index;
    tint = tints(index - 1:index);
end

