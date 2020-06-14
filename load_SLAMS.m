function [SLAMS, database_info] = load_SLAMS(database_name)
    %LOAD_SLAMS Reads SLAMS from file and returns SLAMS structs

    
    SLAMS_db_path = 'C:\Users\carlh\Documents\MATLAB\Exjobb\MMS_SLAMS\SLAMS_database';

    dir_path = [SLAMS_db_path, '\', database_name];

    database_info = read_info(dir_path);

    SLAMS = read_SLAMS(dir_path, database_info);

    if isfile([dir_path, '\search_durations.txt'])
        database_info.search_durations = read_search_durations(dir_path, database_info);
    end

    fclose('all');

    function database_info = read_info(dir_path)
        file_info = fopen([dir_path, '\database_info.txt'], 'r');
        while true
            l = fgetl(file_info);
            if l == -1
                break;
            end

            switch l
                case 'Spacecraft:'
                    database_info.spacecraft = strtrim(fgetl(file_info));
                case 'Time intervals searched:'
                    database_info.search_intervals = read_search_intervals(file_info);
                case 'SLAMS finder settings:'
                    database_info.finder_settings = read_finder_settings(file_info);
                case 'Region class priors:'
                    [classes, priors] = read_class_priors(file_info);
                    database_info.region_classes = classes;
                    database_info.region_class_priors = priors;
                    database_info.n_region_classes = length(classes);
            end
        end
        fclose(file_info);

        function lines = read_indent(fileID)
            lines = {};
            while true
                l = fgetl(fileID);
                if ~strcmp(l(1), sprintf('\t'))
                    fseek(fileID, -(length(l) + 1), 'cof');
                    break;
                else
                    lines = vertcat(lines, strtrim(l)); %#ok<AGROW>
                end
            end
        end
        function tints = read_search_intervals(fileID)
            lines = read_indent(fileID);
            each_date = split(lines, ' - ');
            tints = EpochTT(char(each_date'));
        end
    
        function settings = read_finder_settings(fileID)
            lines = read_indent(fileID);
            settings = split(lines, ' = ');
            settings(:, 2) = cellfun(@eval_or_str, settings(:, 2), 'UniformOutput', false);
            settings = reshape(settings', 1, []);
    
            function x = eval_or_str(x)
                try
                    x = eval(x);
                catch
                    return;
                end
            end
        end
    
        function [classes, priors] = read_class_priors(fileID)
            lines = read_indent(fileID);
            classes_priors = split(lines, ' = ');
            classes = classes_priors(:, 1)';
            priors = cellfun(@str2num, classes_priors(:, 2))';
        end
    end

    function SLAMS = read_SLAMS(dir_path, database_info)
        file_SLAMS = fopen([dir_path, '\SLAMS.csv'], 'r');
        header = fgetl(file_SLAMS);

        formatSpec = '%d %s %s';

        if setting_get(database_info.finder_settings, 'Include_GSE_coords')
            formatSpec = [formatSpec, ' %f %f %f'];
        end

        if setting_get(database_info.finder_settings, 'Include_B_stats')
            formatSpec = [formatSpec, ' %f %f %f %f %f'];
        end

        if setting_get(database_info.finder_settings, 'Include_region_stats')
            n_cols = 2*database_info.n_region_classes + 1;
            formatSpec = [formatSpec, repmat(' %f', [1, n_cols])];

            try
                region_windows = setting_get(database_info.finder_settings, 'Region_time_windows');
            catch
                region_windows = [];
            end

            if ~isempty(region_windows)
                n_windows = length(region_windows);
                n_cols = n_cols*n_windows;
                formatSpec = [formatSpec, repmat(' %f', [1, n_cols])];
            end
        end

        cel = textscan(file_SLAMS, formatSpec, 'Delimiter', ',', 'MultipleDelimsAsOne', true);

        n_SLAMS = length(cel{:, 1});
        cell_converter_tool = repelem(1, n_SLAMS);

        starts = mat2cell(EpochTT(char(cel{:, 2})), cell_converter_tool, 1);
        stops = mat2cell(EpochTT(char(cel{:, 3})), cell_converter_tool, 1);
        ids = mat2cell(cel{:, 1}, cell_converter_tool, 1);
        SLAMS = struct('id', ids, 'start', starts, 'stop', stops);

        if setting_get(database_info.finder_settings, 'Include_GSE_coords')
            col_idx = get_col_idx('x_GSE');
            col_range = col_idx:(col_idx + 2);
            pos = mat2cell([cel{:, col_range}], cell_converter_tool, 3);
            [SLAMS.pos_GSE] = pos{:};
        end

        if setting_get(database_info.finder_settings, 'Include_B_stats')
            col_idx = get_col_idx('B_background_mean');
            B_bg_mean = mat2cell(cel{:, col_idx}, cell_converter_tool, 1);
            B_mean = mat2cell(cel{:, col_idx + 1}, cell_converter_tool, 1);
            B_max = mat2cell(cel{:, col_idx + 2}, cell_converter_tool, 1);
            B_rel_mean = mat2cell(cel{:, col_idx + 3}, cell_converter_tool, 1);
            B_rel_max = mat2cell(cel{:, col_idx + 4}, cell_converter_tool, 1);

            [SLAMS.B_bg_mean] = B_bg_mean{:};
            [SLAMS.B_mean] = B_mean{:};
            [SLAMS.B_max] = B_max{:};
            [SLAMS.B_rel_mean] = B_rel_mean{:};
            [SLAMS.B_rel_max] = B_rel_max{:};           
        end

        if setting_get(database_info.finder_settings, 'Include_region_stats')
            col_idx = get_col_idx(database_info.region_classes{1});
            n_classes = database_info.n_region_classes;

            try
                region_windows = setting_get(database_info.finder_settings, 'Region_time_windows');
                n_windows = length(region_windows);
            catch
                n_windows = 0;
            end

            n_sets = n_windows + 1;
            col_range = col_idx:(col_idx + n_sets*(2*n_classes + 1) - 1);

            tmp = reshape([cel{:, col_range}], n_SLAMS, (2*n_classes + 1), []);
            tmp = permute(tmp, [3, 2, 1]);
            tmp = mat2cell(tmp, n_sets, [n_classes, n_classes, 1], n_SLAMS);
            posteriors = mat2cell(tmp{1}, n_sets, n_classes, cell_converter_tool);
            mahaldists = mat2cell(tmp{2}, n_sets, n_classes, cell_converter_tool);
            logpdfs = mat2cell(tmp{3}, n_sets, 1, cell_converter_tool);

            [SLAMS.region_posterior] = posteriors{:};
            [SLAMS.region_mahaldist] = mahaldists{:};
            [SLAMS.region_logpdf] = logpdfs{:};
        end

        function col_idx = get_col_idx(col_name)
            idx = strfind(header, col_name) - 2;
            idx = idx(1);
            stop_idx = idx;
            while header(idx) ~= '('
                idx = idx - 1;
            end
            start_idx = idx + 1;
            col_idx = str2double(header(start_idx:stop_idx));
        end
    end
end

function map = read_search_durations(dir_path, database_info)
    file_sd = fopen([dir_path, '\search_durations.txt'], 'r');

    formatSpec = '%d,%d = ';

    try
        include_region_stats = setting_get(database_info.finder_settings, 'Include_region_stats');
    catch
        include_region_stats = false;
    end

    if include_region_stats
        wid = database_info.n_region_classes + 1;
        format_add = repmat('%f', [1, wid]);
        format_add = ['[', format_add, ']'];
        formatSpec = [formatSpec, format_add];
    else
        wid = 1;
        format_add = '%f';
        formatSpec = [formatSpec, format_add];
    end
    mat = fscanf(file_sd, formatSpec);
    mat = reshape(mat, wid + 2, [])';

    map = containers.Map;
    [row, col] = size(mat);
    for i = 1:row
        key = sprintf('%d,%d', mat(i, 1), mat(i, 2));
        value = mat(i,3:col);
        map(key) = value;
    end
end