function value = setting_extractor(settings, setting)
%SETTING_EXTRACTOR Summary of this function goes here
%   Detailed explanation goes here
    logic = cellfun(@(x)isequal(x, setting), settings);
    idx = find(logic);
    if ~isempty(idx)
        value = settings{idx + 1};
    else
        error(['Setting ''' setting ''' does not exist.'])
    end
end

