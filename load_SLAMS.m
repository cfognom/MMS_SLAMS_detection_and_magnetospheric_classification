function SLAMS_database = load_SLAMS(database_name)
    %LOAD_SLAMS Reads SLAMS from file and returns SLAMS structs

    
    SLAMS_db_path = 'SLAMS_database';

    dir_path = [SLAMS_db_path, '\', database_name];

    SLAMS_database = read_info(dir_path);

    % read SLAMS
    SLAMS_database.SLAMS_primary = read_SLAMS([dir_path, '\SLAMS_primary.csv'], SLAMS_database);

    if setting_get(SLAMS_database.finder_settings, 'Track_search_durations')
        SLAMS_database.search_durations_primary = read_search_durations([dir_path, '\search_durations_primary.txt'], SLAMS_database);
    end

    if SLAMS_database.include_secondary_intervals
        % read SLAMS where only fgm is active
        original_setting = setting_get(SLAMS_database.finder_settings, 'Include_region_stats');
        SLAMS_database.finder_settings = setting_set(SLAMS_database.finder_settings, 'Include_region_stats', false);

        SLAMS_database.SLAMS_secondary = read_SLAMS([dir_path, '\SLAMS_secondary.csv'], SLAMS_database);
        
        if setting_get(SLAMS_database.finder_settings, 'Track_search_durations')
            SLAMS_database.search_durations_secondary = read_search_durations([dir_path, '\search_durations_secondary.txt'], SLAMS_database);
        end

        SLAMS_database.finder_settings = setting_set(SLAMS_database.finder_settings, 'Include_region_stats', original_setting);
    end

    fclose('all');

    function SLAMS_database = read_info(dir_path)
        file_info = fopen([dir_path, '\database_info.txt'], 'r');
        while true
            l = fgetl(file_info);
            if l == -1
                break;
            end

            switch l
                case 'Spacecraft:'
                    SLAMS_database.spacecraft = strtrim(fgetl(file_info));
                case 'Time intervals searched:'
                    SLAMS_database.search_intervals = read_search_intervals(file_info);
                case 'SLAMS finder settings:'
                    SLAMS_database.finder_settings = read_finder_settings(file_info);
                case 'Region class priors:'
                    [classes, priors] = read_class_priors(file_info);
                    SLAMS_database.region_classes = classes;
                    SLAMS_database.region_class_priors = priors;
                    SLAMS_database.n_region_classes = length(classes);
                case 'Other:'
                    SLAMS_database.include_secondary_intervals = read_other(file_info);
            end
        end
        fclose(file_info);
        
        function lines = read_indent(fileID)
            lines = {};
            while true
                l = fgetl(fileID);
                if ~strcmp(l(1), sprintf('\t'))
                    fseek(fileID, -(length(l) + 2), 'cof');
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

        function include_secondary_intervals = read_other(fileID)
            lines = read_indent(fileID);
            setting_value = split(lines, ' = ')';
            if strcmp(setting_value{1,1}, 'Include_secondary_intervals')
                include_secondary_intervals = eval(setting_value{1,2});
            end
        end
    end

    function SLAMS = read_SLAMS(SLAMS_file, SLAMS_database)
        file_SLAMS = fopen(SLAMS_file, 'r');
        header = fgetl(file_SLAMS);

        formatSpec = '%d %s %s';

        if setting_get(SLAMS_database.finder_settings, 'Include_GSE_coords')
            formatSpec = [formatSpec, ' %f %f %f'];
        end

        if setting_get(SLAMS_database.finder_settings, 'Include_B_stats')
            formatSpec = [formatSpec, ' %f %f %f %f %f'];
        end

        if setting_get(SLAMS_database.finder_settings, 'Include_region_stats')
            n_cols = 2*SLAMS_database.n_region_classes + 1;
            formatSpec = [formatSpec, repmat(' %f', [1, n_cols])];

            try
                region_windows = setting_get(SLAMS_database.finder_settings, 'Region_time_windows');
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

        if setting_get(SLAMS_database.finder_settings, 'Include_GSE_coords')
            col_idx = get_col_idx('x_GSE');
            col_range = col_idx:(col_idx + 2);
            pos = mat2cell([cel{:, col_range}], cell_converter_tool, 3);
            [SLAMS.pos_GSE] = pos{:};
        end

        if setting_get(SLAMS_database.finder_settings, 'Include_B_stats')
            col_idx = get_col_idx('B0_mean');
            B0_mean = mat2cell(cel{:, col_idx}, cell_converter_tool, 1);
            B_mean = mat2cell(cel{:, col_idx + 1}, cell_converter_tool, 1);
            B_max = mat2cell(cel{:, col_idx + 2}, cell_converter_tool, 1);
            B_rel_mean = mat2cell(cel{:, col_idx + 3}, cell_converter_tool, 1);
            B_rel_max = mat2cell(cel{:, col_idx + 4}, cell_converter_tool, 1);

            [SLAMS.B0_mean] = B0_mean{:};
            [SLAMS.B_mean] = B_mean{:};
            [SLAMS.B_max] = B_max{:};
            [SLAMS.B_rel_mean] = B_rel_mean{:};
            [SLAMS.B_rel_max] = B_rel_max{:};           
        end

        if setting_get(SLAMS_database.finder_settings, 'Include_region_stats')
            col_idx = get_col_idx(SLAMS_database.region_classes{1});
            n_classes = SLAMS_database.n_region_classes;

            try
                region_windows = setting_get(SLAMS_database.finder_settings, 'Region_time_windows');
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

function map = read_search_durations(file_name, SLAMS_database)
    file_sd = fopen(file_name, 'r');

    formatSpec = '%d,%d = ';

    try
        include_region_stats = setting_get(SLAMS_database.finder_settings, 'Include_region_stats');
    catch
        include_region_stats = false;
    end

    if include_region_stats
        wid = SLAMS_database.n_region_classes + 1;
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