sc = 'MMS1';

plot_prediction = false;

save_file_name = 'identified_SLAMS.csv';

tints_search = search_intervals();

settings = {
    'Include_B_stats', true, ...
    'Include_region_stats', true, ...
    'Region_time_windows', [15, 30, 60, 2*60, 4*60, 8*60], ...
    'Include_GSE_coords', true, ...
    'Extra_load_time', 4*60, ...
    'SLAMS_B_bg_method', 'median', ...
    'SLAMS_B_bg_window', 30, ...
    'SLAMS_threshold', 2, ...
    'SLAMS_min_duration', 0
};

main(sc, tints_search, settings, save_file_name, plot_prediction);

function main(sc, tints_search, settings, save_file_name, plot_prediction)
    finder = SLAMS_finder('Spacecraft', sc, 'Show_train_progress', false);
    
    tints_active = get_tints_active(tints_search);
    
    tints_valid = remove_short_tints(tints_active, 2*setting_get(settings, 'Extra_load_time'));
    
    n_tints = length(tints_valid)/2;
    
    fileID = fopen(save_file_name, 'w');
    write_info(fileID, settings, tints_search, finder);
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
        current_id = write_SLAMS(fileID, SLAMS, settings, current_id, finder);
    end
    disp('Done!')
    % fclose(fileID);
    fclose('all');
end

function write_info(fileID, settings, tints_search, finder)

    sc_str = {'Spacecraft:', finder.sc};
    sc_str = [strjoin(sc_str, '\n\t'), '\n'];
    fprintf(fileID, sc_str);

    n_t_seach_user = length(tints_search);
    search_tint_str = mat2cell(tints_search.utc, repelem(1, n_t_seach_user), 30);
    search_tint_str = reshape(search_tint_str, 2, [])';
    search_tint_str = join(search_tint_str, ' - ');
    search_tint_str = vertcat({'Time intervals searched:'}, search_tint_str);
    search_tint_str = [strjoin(search_tint_str, '\n\t'), '\n'];
    fprintf(fileID, search_tint_str);

    finder_settings_str = construct_settings_string(settings);
    finder_settings_str = vertcat({'SLAMS finder settings:'}, finder_settings_str);
    finder_settings_str = [strjoin(finder_settings_str, '\n\t'), '\n'];
    fprintf(fileID, finder_settings_str);

    if setting_get(settings, 'Include_region_stats')
        priors_str = join([finder.classes', split(num2str(finder.priors))], ' = ');
        priors_str = vertcat({'Region class priors:'}, priors_str);
        priors_str = [strjoin(priors_str, '\n\t'), '\n'];
        fprintf(fileID, priors_str);
    end

    fprintf(fileID, 'Identified SLAMS:\n');

    header_str = construct_header(settings);
    fprintf(fileID, header_str);
end

function current_id = write_SLAMS(fileID, SLAMS, settings, current_id, finder)
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
            cell_converter_tool2 = repelem(1, finder.n_classes);
            n_cols = (2*finder.n_classes + 1);
            add1 = mat2cell(vertcat(SLAMS.region_posterior), cell_converter_tool, cell_converter_tool2);
            add2 = mat2cell(vertcat(SLAMS.region_mahaldist), cell_converter_tool, cell_converter_tool2);
            add3 = mat2cell(vertcat(SLAMS.region_logpdf), cell_converter_tool, 1);
            cel = horzcat(cel, add1, add2, add3);
            str = [str, repmat(', %f', [1, n_cols])];
            
            windows = setting_get(settings, 'Region_time_windows');
            if ~isempty(windows)
                n_windows = length(windows);
                n_cols = n_cols*n_windows;
                cell_converter_tool3 = repelem(1, n_cols);
                
                posteriors = cat(3, SLAMS.region_posterior_windows);
                posteriors = permute(posteriors, [3, 2, 1]);
                mahaldists = cat(3, SLAMS.region_mahaldist_windows);
                mahaldists = permute(mahaldists, [3, 2, 1]);
                logpdfs = [SLAMS.region_logpdf_windows];
                logpdfs = permute(logpdfs, [2, 3, 1]);
                tmp = cat(2, posteriors, mahaldists, logpdfs);
                tmp = reshape(tmp, n_SLAMS, []);

                add1 = mat2cell(tmp, cell_converter_tool, cell_converter_tool3);

                cel = horzcat(cel, add1);
                str = [str, repmat(', %f', [1, n_cols])];
            end
        end

        cel = cel';
        str = [str, '\n'];
        fprintf(fileID, str, cel{:});
    end
end

function header = construct_header(settings)
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
        header_add = {
            'MSP_posterior', 'SW_posterior', 'MSH_posterior', 'FS_posterior', ...
            'MSH_mahaldist', 'SW_mahaldist', 'MSP_mahaldist', 'FS_posterior', ...
            'logpdf'
        };
        header = horzcat(header, header_add);

        windows = setting_get(settings, 'Region_time_windows');
        if ~isempty(windows)
            str_windows = split(num2str(windows))';
            for w = str_windows
                for t = {'posterior', 'mahaldist'}
                    for r = {'MSP', 'SW', 'MSH', 'FS'}
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