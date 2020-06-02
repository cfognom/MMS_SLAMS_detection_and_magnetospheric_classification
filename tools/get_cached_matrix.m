function matrix = get_cached_matrix(name, func)
%GET_CACHED_MATRIX Summary of this function goes here
%   Detailed explanation goes here

    file_path = ['cache/', name, '.txt'];
    if isfile(file_path)
        fprintf('Loading existing file ''%s''.\n', file_path)
        matrix = readmatrix(file_path);
        % tints = EpochTT(int64(raw));
    else
        fprintf('Generating new file ''%s''...\n', file_path)
        matrix = func();
        writematrix(matrix, file_path);
    end
end

