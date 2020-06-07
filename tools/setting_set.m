function settings = setting_set(settings, setting, value)
    %SETTING_SET Summary of this function goes here
    %   Detailed explanation goes here
        logic = cellfun(@(x)isequal(x, setting), settings);
        idx = find(logic);
        if ~isempty(idx)
            settings{idx + 1} = value;
        else
            error(['Setting ''' setting ''' does not exist.'])
        end
    end

