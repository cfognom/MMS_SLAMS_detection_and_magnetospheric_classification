function settings = setting_setter(settings, setting, value)
    %SETTING_SETTER Summary of this function goes here
    %   Detailed explanation goes here
        logic = cellfun(@(x)isequal(x, setting), settings);
        idx = find(logic);
        if ~isempty(idx)
            settings{idx + 1} = value;
        else
            error(['Setting ''' setting ''' does not exist.'])
        end
    end

