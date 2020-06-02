function splits = split_array(array, n)
%SPLIT_ARRAY Summary of this function goes here
%   Detailed explanation goes here
    l = length(array);
    sub_l = l/n;
    splits = cell(1, n);
    for i = 1:n
        splits{i} = array(int32((i-1)*sub_l) + 1:int32(i*sub_l));
    end
end

