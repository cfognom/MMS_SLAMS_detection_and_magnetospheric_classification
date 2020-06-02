function tints = union_tints(varargin)
%UNION_TINTS Summary of this function goes here
%   Detailed explanation goes here
    % n_vars = length(varargin);
    starts = cellfun(@(x) x(1:2:end), varargin, 'UniformOutput', false);
    starts = [starts{:}];
    stops = cellfun(@(x) x(2:2:end), varargin, 'UniformOutput', false);
    stops = [stops{:}];
    type_vector = [ones(1, length(starts)), -ones(1, length(stops))];
    all = [starts, stops];
    [all, sort_order] = sort(all);
    type_vector = [0, type_vector(sort_order)];
    tints = all(diff(cumsum(type_vector) >= 1) ~= 0);
end

