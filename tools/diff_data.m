function dif = diff_data(data, dim)
%DIFFERENTIATE_DATA Summary of this function goes here
%   Detailed explanation goes here
    if dim == 2
        data = data';
    end
    dif = diff(data);
    dif1 = [dif; dif(end, :)];
    dif2 = [dif(1, :); dif];
    dif = (dif1 + dif2)/2;
    if dim == 2
        dif = dif';
    end
end

