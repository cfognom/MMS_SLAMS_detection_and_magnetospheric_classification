function out = get_cached_var(var_name)
%GET_CACHED_VAR Summary of this function goes here
%   Detailed explanation goes here
    file_path = ['cache/', var_name];
    if isfile(file_path)
        out = load(file_path);
    else
        out = [];
    %     if nargin > 1
    %         out = func_generate();
    %         save(file_path, var_name)
        % else
            % error('Could not find variable ''%s'' in cache.', var_name);
    %     end
    end
end

