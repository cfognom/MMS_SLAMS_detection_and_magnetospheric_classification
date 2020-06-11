function [SLAMS, dataset_info] = load_SLAMS(file, varargin)
    %LOAD_SLAMS Reads SLAMS from file and returns SLAMS structs
    % Can filter SLAMS using a the following name value pairs:
    %   'ID', [Numeric array of ids]
    %   'After_date', [UTC string]
    %   'Before_date', [UTC string]
    %   'Min_duration', [Numeric value in seconds]
    %   'Max_duration', [Numeric value in seconds]
    %   'GSE_filter', [Handle to a function with signature:
    %       bool = myFunc(x_GSE, y_GSE, z_GSE).
    %       If bool is true, SLAMS is included]
    %   'B_filter', [Handle to a function with signature:
    %       bool = myFunc(B_background_mean, B_mean, B_max, B_relative_mean, B_relative_max).
    %       If bool is true, SLAMS is included]
    %   'Region_filter', [Handle to a function with signature:
    %       bool = myFunc(B_background_mean, B_mean, B_max, B_relative_mean, B_relative_max).
    %       If bool is true, SLAMS is included]
    
    
    % load_filter = {
    %     'ID', []
    %     'After_date', []
    %     'Before_date', []
    %     'Min_duration', []
    %     'Max_duration', []
    %     'GSE_filter', []
    %     'B_filter', []
    %     'Region_filter', []
    % };

    p = inputParser;
    addRequired(p, 'file');
    addParameter(p, 'ID', [])
    addParameter(p, 'After_date', [])
    addParameter(p, 'Before_date', [])
    addParameter(p, 'Min_duration', [])
    addParameter(p, 'Max_duration', [])
    addParameter(p, 'GSE_filter', [])
    addParameter(p, 'B_filter', [])
    addParameter(p, 'Region_filter', [])
    parse(p, file, varargin{:})
    r = p.Results;

    fileID = fopen(file, 'r');
    
    dataset_info = read_info(fileID);
    
    SLAMS = read_SLAMS(fileID, dataset_info);
    
    % if isempty(r.ID)

    % else

    % end

    fclose('all');

    function dataset_info = read_info(fileID)
        while true
            l = fgetl(fileID);

            switch l
                case 'Spacecraft:'
                    dataset_info.spacecraft = strtrim(fgetl(fileID));
                case 'Time intervals searched:'
                    dataset_info.search_intervals = read_search_intervals(fileID);
                case 'SLAMS finder settings:'
                    dataset_info.finder_settings = read_finder_settings(fileID);
                case 'Region class priors:'
                    [classes, priors] = read_class_priors(fileID);
                    dataset_info.region_classes = classes;
                    dataset_info.region_class_priors = priors;
                    dataset_info.n_region_classes = length(classes);
                case 'Identified SLAMS:'
                    return;
            end
        end

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

    function SLAMS = read_SLAMS(fileID, dataset_info)
        header = fgetl(fileID);

        formatSpec = '%d %s %s';

        if setting_get(dataset_info.finder_settings, 'Include_GSE_coords')
            formatSpec = [formatSpec, ' %f %f %f'];
        end

        if setting_get(dataset_info.finder_settings, 'Include_B_stats')
            formatSpec = [formatSpec, ' %f %f %f %f %f'];
        end

        if setting_get(dataset_info.finder_settings, 'Include_region_stats')
            n_cols = 2*dataset_info.n_region_classes + 1;
            formatSpec = [formatSpec, repmat(' %f', [1, n_cols])];

            try
                region_windows = setting_get(dataset_info.finder_settings, 'Region_time_windows');
            catch
                region_windows = [];
            end

            if ~isempty(region_windows)
                n_windows = length(region_windows);
                n_cols = n_cols*n_windows;
                formatSpec = [formatSpec, repmat(' %f', [1, n_cols])];
            end
        end

        cel = textscan(fileID, formatSpec, 'Delimiter', ',', 'MultipleDelimsAsOne', true);

        n_SLAMS = length(cel{:, 1});
        cell_converter_tool = repelem(1, n_SLAMS);

        starts = mat2cell(EpochTT(char(cel{:, 2})), cell_converter_tool, 1);
        stops = mat2cell(EpochTT(char(cel{:, 3})), cell_converter_tool, 1);
        ids = mat2cell(cel{:, 1}, cell_converter_tool, 1);
        SLAMS = struct('id', ids, 'start', starts, 'stop', stops);

        if setting_get(dataset_info.finder_settings, 'Include_GSE_coords')
            col_idx = get_col_idx('x_GSE');
            col_range = col_idx:(col_idx + 2);
            pos = mat2cell([cel{:, col_range}], cell_converter_tool, 3);
            [SLAMS.pos_GSE] = pos{:};
        end

        if setting_get(dataset_info.finder_settings, 'Include_B_stats')
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

        if setting_get(dataset_info.finder_settings, 'Include_region_stats')
            col_idx = get_col_idx(dataset_info.region_classes{1});
            n_classes = dataset_info.n_region_classes;
            col_range_posteriors = col_idx:(col_idx + n_classes - 1);
            col_range_mahaldist = col_range_posteriors + n_classes;
            col_logpdf = col_range_mahaldist(end) + 1;

            posteriors = mat2cell([cel{:, col_range_posteriors}], cell_converter_tool, n_classes);
            mahaldist = mat2cell([cel{:, col_range_mahaldist}], cell_converter_tool, n_classes);
            logpdf = mat2cell([cel{:, col_logpdf}], cell_converter_tool, 1);

            [SLAMS.region_posterior] = posteriors{:};
            [SLAMS.region_mahaldist] = mahaldist{:};
            [SLAMS.region_logpdf] = logpdf{:};

            try
                region_windows = setting_get(dataset_info.finder_settings, 'Region_time_windows');
            catch
                region_windows = [];
            end

            if ~isempty(region_windows)
                n_windows = length(region_windows);
                col_idx = col_logpdf + 1;
                col_range = col_idx:(col_idx + n_windows*(2*n_classes + 1) - 1);

                tmp = reshape([cel{:, col_range}], n_SLAMS, (2*n_classes + 1), []);
                tmp = permute(tmp, [3, 2, 1]);
                tmp = mat2cell(tmp, n_windows, [n_classes, n_classes, 1], n_SLAMS);
                posterior_windows = mat2cell(tmp{1}, n_windows, n_classes, cell_converter_tool);
                mahaldist_windows = mat2cell(tmp{2}, n_windows, n_classes, cell_converter_tool);
                logpdf_windows = mat2cell(tmp{3}, n_windows, 1, cell_converter_tool);

                [SLAMS.region_posterior_windows] = posterior_windows{:};
                [SLAMS.region_mahaldist_windows] = mahaldist_windows{:};
                [SLAMS.region_logpdf_windows] = logpdf_windows{:};
            end
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
