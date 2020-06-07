function value = setting_get(settings, setting)
%SETTING_GET Summary of this function goes here
%   Detailed explanation goes here
    logic = cellfun(@(x)isequal(x, setting), settings);
    idx = find(logic);
    if ~isempty(idx)
        value = settings{idx + 1};
    else
        error(['Setting ''' setting ''' does not exist.'])
    end
end

