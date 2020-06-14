sc = 'MMS1';

plot_prediction = false;

SLAMS_db_path = 'C:\Users\carlh\Documents\MATLAB\Exjobb\MMS_SLAMS\SLAMS_database';

database_name = 'identified_SLAMS';

tints_search = search_intervals();

settings = {
    'Include_B_stats', true, ...
    'Include_region_stats', true, ...
    'Region_time_windows', [15, 30, 60, 2*60, 4*60, 8*60], ...
    'Include_GSE_coords', true, ...
    'Track_search_durations', true, ...
    'Extra_load_time', 4*60, ...
    'SLAMS_B_bg_method', 'median', ...
    'SLAMS_B_bg_window', 60, ...
    'SLAMS_threshold', 2, ...
    'SLAMS_min_duration', 0
};

main(sc, tints_search, settings, SLAMS_db_path, database_name, plot_prediction);

function main(sc, tints_search, settings, SLAMS_db_path, database_name, plot_prediction)
    dir_path = [SLAMS_db_path, '\', database_name];
    if ~exist(dir_path, 'dir')
        mkdir(dir_path);
    end

    finder = SLAMS_finder('Spacecraft', sc, 'Show_region_classifier_steps', false);
    
    tints_active = get_tints_active(tints_search);
    
    tints_valid = remove_short_tints(tints_active, 2*setting_get(settings, 'Extra_load_time'));
    
    n_tints = length(tints_valid)/2;
    
    write_info(dir_path, settings, tints_search, finder);
    file_SLAMS = fopen([dir_path, '\SLAMS.csv'], 'w');
    header_str = construct_header(settings, finder);
    fprintf(file_SLAMS, header_str);
    current_id = 1;
    for i = 1:n_tints
        fprintf('Looking for SLAMS in interval %u/%u\n', i, n_tints);
        tint = select_tint(tints_valid, i);
        extra_load_time = setting_get(settings, 'Extra_load_time');
        tint = [tint(1) + extra_load_time, tint(2) + -extra_load_time];
        SLAMS = finder.evaluate(tint, settings{:});
        if plot_prediction
            finder.plot_prediction(tint);
        end
        current_id = write_SLAMS(file_SLAMS, SLAMS, settings, current_id, finder);
        % break;
    end

    try
        track_search_durations = setting_get(settings, 'Track_search_durations');
    catch
        track_search_durations = false;
    end
    
    if track_search_durations
        write_search_durations(dir_path, finder);
    end

    disp('Done!')
    % fclose(file_SLAMS);
    fclose('all');
end

function write_info(dir_path, settings, tints_search, finder)
    file_info = fopen([dir_path, '\database_info.txt'], 'w');

    sc_str = {'Spacecraft:', finder.sc};
    sc_str = [strjoin(sc_str, '\n\t'), '\n'];
    fprintf(file_info, sc_str);

    n_t_seach_user = length(tints_search);
    search_tint_str = mat2cell(tints_search.utc, repelem(1, n_t_seach_user), 30);
    search_tint_str = reshape(search_tint_str, 2, [])';
    search_tint_str = join(search_tint_str, ' - ');
    search_tint_str = vertcat({'Time intervals searched:'}, search_tint_str);
    search_tint_str = [strjoin(search_tint_str, '\n\t'), '\n'];
    fprintf(file_info, search_tint_str);

    finder_settings_str = construct_settings_string(settings);
    finder_settings_str = vertcat({'SLAMS finder settings:'}, finder_settings_str);
    finder_settings_str = [strjoin(finder_settings_str, '\n\t'), '\n'];
    fprintf(file_info, finder_settings_str);

    if setting_get(settings, 'Include_region_stats')
        priors_str = join([finder.classes', split(num2str(finder.priors))], ' = ');
        priors_str = vertcat({'Region class priors:'}, priors_str);
        priors_str = [strjoin(priors_str, '\n\t'), '\n'];
        fprintf(file_info, priors_str);
    end

    fclose(file_info);
end

function current_id = write_SLAMS(file_SLAMS, SLAMS, settings, current_id, finder)
    if ~isempty(SLAMS)
        n_SLAMS = length(SLAMS);
        fprintf('Found %u SLAMS!\n', n_SLAMS);
        id = int64((current_id:(current_id + n_SLAMS - 1))');
        current_id = current_id + n_SLAMS;
        
        cell_converter_tool = repelem(1, n_SLAMS);

        cel = mat2cell(id, cell_converter_tool, 1);
        str = '%u';

        SLAMS_starts = [SLAMS.start];
        SLAMS_stops = [SLAMS.stop];
        add1 = mat2cell(SLAMS_starts.utc, cell_converter_tool, 30);
        add2 = mat2cell(SLAMS_stops.utc, cell_converter_tool, 30);
        cel = horzcat(cel, add1, add2);
        str = [str, ', %s, %s'];

        if setting_get(settings, 'Include_GSE_coords')
            add1 = mat2cell(vertcat(SLAMS.pos_GSE), cell_converter_tool, [1, 1, 1]);
            cel = horzcat(cel, add1);
            str = [str, ', %f, %f, %f'];
        end
        
        if setting_get(settings, 'Include_B_stats')
            add1 = mat2cell(vertcat(SLAMS.B_bg_mean), cell_converter_tool, 1);
            add2 = mat2cell(vertcat(SLAMS.B_mean), cell_converter_tool, 1);
            add3 = mat2cell(vertcat(SLAMS.B_max), cell_converter_tool, 1);
            add4 = mat2cell(vertcat(SLAMS.B_rel_mean), cell_converter_tool, 1);
            add5 = mat2cell(vertcat(SLAMS.B_rel_max), cell_converter_tool, 1);
            cel = horzcat(cel, add1, add2, add3, add4, add5);
            str = [str, ', %f, %f, %f, %f, %f'];
        end
        
        if setting_get(settings, 'Include_region_stats')
            
            try
                region_windows = setting_get(settings, 'Region_time_windows');
                n_windows = length(region_windows);
            catch
                n_windows = 0;
            end
            
            n_sets = n_windows + 1;
            n_cols = n_sets*(2*finder.n_classes + 1);
            cell_converter_tool2 = repelem(1, n_cols);

            posteriors = permute(cat(3, SLAMS.region_posterior), [3, 2, 1]);
            mahaldists = permute(cat(3, SLAMS.region_mahaldist), [3, 2, 1]);
            logpdfs = permute(cat(3, SLAMS.region_logpdf), [3, 2, 1]);
            tmp = horzcat(posteriors, mahaldists, logpdfs);
            tmp = reshape(tmp, n_SLAMS, []);

            add1 = mat2cell(tmp, cell_converter_tool, cell_converter_tool2);
            cel = horzcat(cel, add1);
            str = [str, repmat(', %f', [1, n_cols])];
        end

        cel = cel';
        str = [str, '\n'];
        fprintf(file_SLAMS, str, cel{:});
    end
end

function write_search_durations(dir_path, finder)
    file_sd = fopen([dir_path, '\search_durations.txt'], 'w');

    map = finder.search_durations;
    % n_bins = map.Count;
    k = keys(map);
    for i = k
        fprintf(file_sd, [i{:}, ' = ', mat2str(map(i{:})), '\n']);
    end
    fclose(file_sd);
end

function header = construct_header(settings, finder)
    header = {'id', 'start_UTC', 'stop_UTC'};

    if setting_get(settings, 'Include_GSE_coords')
        header_add = {'x_GSE', 'y_GSE', 'z_GSE'};
        header = horzcat(header, header_add);
    end

    if setting_get(settings, 'Include_B_stats')
        header_add = {'B_background_mean', 'B_mean', 'B_max', 'B_relative_mean', 'B_relative_max'};
        header = horzcat(header, header_add);
    end

    if setting_get(settings, 'Include_region_stats')
        for t = {'posterior', 'mahaldist'}
            for r = finder.classes
                header = horzcat(header, {[r{:} '_' t{:}]}); %#ok<AGROW>
            end
        end
        header = horzcat(header, {'logpdf'});

        windows = setting_get(settings, 'Region_time_windows');
        if ~isempty(windows)
            str_windows = split(num2str(windows))';
            for w = str_windows
                for t = {'posterior', 'mahaldist'}
                    for r = finder.classes
                        header = horzcat(header, {[r{:} '_' t{:} '_' w{:} 'sec_mean']}); %#ok<AGROW>
                    end
                end
                header = horzcat(header, {['logpdf_' w{:} 'sec_mean']}); %#ok<AGROW>
            end
        end
    end

    n_cols = length(header);
    column_id = split(num2str(1:n_cols));
    parenthesis = join([repelem({'('}, n_cols)', column_id, repelem({')'}, n_cols)'], '');
    header = join(join([parenthesis, header'], ''), ', ');
    header = header{:};
    header = [header, '\n'];
end