function tints = intersect_tints(varargin)
%INTERSECT_INT Summary of this function goes here
%   Detailed explanation goes here

    % in = horzcat(varargin{:});
    % n = length(varargin);
    % starts = [in.start];
    % stops = [in.stop];
    % type_vector = [ones(1, length(starts)), -ones(1, length(stops))];
    % all = [starts stops];
    % [all, id] = sort(all, 'ascend');
    % type_vector = type_vector(id);
    % counter = 0;
    % out_counter = 1;
    % for i = 1:length(type_vector)
    %     counter = counter + type_vector(i);
    %     if counter == n && type_vector(i) == 1
    %         saved_date = all(i);
    %     elseif counter == n - 1 && type_vector(i) == - 1
    %         out(out_counter) = datetime_int(saved_date, all(i));
    %         out_counter = out_counter + 1;
    %     end
    % end

    n_vars = length(varargin);
    starts = cellfun(@(x) x(1:2:end), varargin, 'UniformOutput', false);
    starts = [starts{:}];
    stops = cellfun(@(x) x(2:2:end), varargin, 'UniformOutput', false);
    stops = [stops{:}];
    type_vector = [ones(1, length(starts)), -ones(1, length(stops))];
    all = [starts, stops];
    [all, sort_order] = sort(all);
    type_vector = [0, type_vector(sort_order)];
    tints = all(diff(cumsum(type_vector) >= n_vars) ~= 0);
end

