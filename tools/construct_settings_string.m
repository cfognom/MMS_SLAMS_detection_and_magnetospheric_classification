function str = construct_settings_string(settings)
%CONSTRUCT_SETTINGS_STRING Summary of this function goes here
%   Detailed explanation goes here
    l = numel(settings);
    str = cell(2, l/2);
    for i = 1:l
        elem = settings{i};
        if ischar(elem)
            str{i} = elem;
        else 
            str{i} = mat2str(elem);
        end
    end
    str = join(str', ' = ');
    % str = strjoin(str, '\n');
end

