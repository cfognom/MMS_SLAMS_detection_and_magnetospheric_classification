%% Start
db_path = 'C:\Users\carlh\Documents\MATLAB\Exjobb\Carl\data\';
cache_path = 'cache\';
sc = 'mms1';

tints_search_user = get_tints_user();

tints_active = get_tints_active(tints_search_user);

finder = SLAMS_finder('Show_train_progress', true);

settings = {
    'Include_B_stats', true, ...
    'Include_region_stats', true, ...
    'Region_time_windows', [15, 30, 60, 2*60, 4*60, 8*60], ...
    'Include_GSE_coords', true, ...
    'Extra_load_time', 4*60, ...
    'SLAMS_B_bg_method', 'median', ...
    'SLAMS_B_bg_window', 60, ...
    'SLAMS_threshold', 2
};

n_tints = length(tints_active)/2;

fileID = fopen('identified_SLAMS.csv', 'w');
print_info(fileID, settings, tints_search_user);
% plot_prediction = true;
plot_prediction = false;
current_id = 1;
for i = 23:n_tints
% for i = 9:n_tints
% for i = 700:n_tints
    fprintf('Looking for SLAMS in interval %u/%u\n', i, n_tints);
    tint = select_tint(tints_active, i);
    extra_load_time = setting_extractor(settings, 'Extra_load_time');
    tint = [tint(1) + extra_load_time, tint(2) + -extra_load_time];
    SLAMS = finder.evaluate(tint, settings{:});
    if plot_prediction
        finder.plot_prediction(tint);
    end
    current_id = print_SLAMS(fileID, SLAMS, settings, current_id);
    % break;
end
disp('Done!')
fclose(fileID);

function print_info(fileID, settings, tints_search_user)

    n_t_seach_user = length(tints_search_user);
    search_tint_str = mat2cell(tints_search_user.utc, repelem(1, n_t_seach_user), 30);
    search_tint_str = reshape(search_tint_str, 2, [])';
    search_tint_str = join(search_tint_str, ' - ');
    search_tint_str = vertcat({'Time intervals searched:'}, search_tint_str);
    search_tint_str = [strjoin(search_tint_str, '\n\t'), '\n'];
    fprintf(fileID, search_tint_str);

    settings_str = construct_settings_string(settings);
    settings_str = vertcat({'Settings:'}, settings_str);
    settings_str = [strjoin(settings_str, '\n\t'), '\n'];
    fprintf(fileID, settings_str);

    fprintf(fileID, 'Identified SLAMS:\n');

    header_str = construct_header(settings);
    fprintf(fileID, header_str);
end

function current_id = print_SLAMS(fileID, SLAMS, settings, current_id)
    if ~isempty(SLAMS)
        n_SLAMS = length(SLAMS);
        fprintf('Found %u SLAMS!\n', n_SLAMS);
        id = int64((current_id:(current_id + n_SLAMS - 1))');
        current_id = current_id + n_SLAMS;
        
        % cel = cell(n_SLAMS, n_cols);
        cell_converter_tool = repelem(1, n_SLAMS);

        cel = mat2cell(id, cell_converter_tool, 1);
        str = '%u';

        SLAMS_starts = [SLAMS.start];
        SLAMS_stops = [SLAMS.stop];
        add1 = mat2cell(SLAMS_starts.utc, cell_converter_tool, 30);
        add2 = mat2cell(SLAMS_stops.utc, cell_converter_tool, 30);
        add3 = mat2cell(SLAMS_starts.epoch, cell_converter_tool, 1);
        add4 = mat2cell(SLAMS_stops.epoch, cell_converter_tool, 1);
        cel = horzcat(cel, add1, add2, add3, add4);
        str = [str, ', %s, %s, %u, %u'];

        if setting_extractor(settings, 'Include_GSE_coords')
            add1 = mat2cell(vertcat(SLAMS.pos_GSE), cell_converter_tool, [1, 1, 1]);
            cel = horzcat(cel, add1);
            str = [str, ', %f, %f, %f'];
        end
        
        if setting_extractor(settings, 'Include_B_stats')
            add1 = mat2cell(vertcat(SLAMS.B_bg_mean), cell_converter_tool, 1);
            add2 = mat2cell(vertcat(SLAMS.B_mean), cell_converter_tool, 1);
            add3 = mat2cell(vertcat(SLAMS.B_max), cell_converter_tool, 1);
            add4 = mat2cell(vertcat(SLAMS.B_rel_mean), cell_converter_tool, 1);
            add5 = mat2cell(vertcat(SLAMS.B_rel_max), cell_converter_tool, 1);
            cel = horzcat(cel, add1, add2, add3, add4, add5);
            str = [str, ', %f, %f, %f, %f, %f'];
        end
        
        if setting_extractor(settings, 'Include_region_stats')
            add1 = mat2cell(vertcat(SLAMS.region_posterior), cell_converter_tool, [1, 1, 1]);
            add2 = mat2cell(vertcat(SLAMS.region_mahaldist), cell_converter_tool, [1, 1, 1]);
            add3 = mat2cell(vertcat(SLAMS.region_logpdf), cell_converter_tool, 1);
            cel = horzcat(cel, add1, add2, add3);
            str = [str, ', %f, %f, %f, %f, %f, %f, %f'];
            
            windows = setting_extractor(settings, 'Region_time_windows');
            if ~isempty(windows)
                n_windows = length(windows);
                cell_converter_tool2 = repelem(1, 7*n_windows);
                
                posteriors = cat(3, SLAMS.region_posterior_windows);
                posteriors = permute(posteriors, [3, 2, 1]);
                mahaldists = cat(3, SLAMS.region_mahaldist_windows);
                mahaldists = permute(mahaldists, [3, 2, 1]);
                logpdfs = [SLAMS.region_logpdf_windows];
                logpdfs = permute(logpdfs, [2, 3, 1]);
                tmp = cat(2, posteriors, mahaldists, logpdfs);
                tmp = reshape(tmp, n_SLAMS, []);

                add1 = mat2cell(tmp, cell_converter_tool, cell_converter_tool2);

                cel = horzcat(cel, add1);
                str = [str, repmat(', %f', [1, 7*n_windows])];
            end
        end

        cel = cel';
        str = [str, '\n'];
        fprintf(fileID, str, cel{:});
    end
end

function header = construct_header(settings)
    header = {'id', 'start_UTC', 'stop_UTC', 'start_EpochTT', 'stop_EpochTT'};

    if setting_extractor(settings, 'Include_GSE_coords')
        header_add = {'x_GSE', 'y_GSE', 'z_GSE'};
        header = horzcat(header, header_add);
    end

    if setting_extractor(settings, 'Include_B_stats')
        header_add = {'B_background_mean', 'B_mean', 'B_max', 'B_relative_mean', 'B_relative_max'};
        header = horzcat(header, header_add);
    end

    if setting_extractor(settings, 'Include_region_stats')
        header_add = {
            'MSH_posterior', 'SW_posterior', 'MSP_posterior', ...
            'MSH_mahaldist', 'SW_mahaldist', 'MSP_mahaldist', ...
            'logpdf'
        };
        header = horzcat(header, header_add);

        windows = setting_extractor(settings, 'Region_time_windows');
        if ~isempty(windows)
            str_windows = split(num2str(windows))';
            for w = str_windows
                for t = {'posterior', 'mahaldist'}
                    for r = {'MSH', 'SW', 'MSP'}
                        header = horzcat(header, {[r{:} '_' t{:} '_' w{:} 'sec_mean']}); %#ok<AGROW>
                    end
                end
                header = horzcat(header, {['logpdf_' w{:} 'sec_mean']}); %#ok<AGROW>
            end
        end
    end

        % ...
        % 'MSH_posterior_15sec_mean', 'SW_posterior_15sec_mean', 'MSP_posterior_15sec_mean', ...
        % 'MSH_posterior_30sec_mean', 'SW_posterior_30sec_mean', 'MSP_posterior_30sec_mean', ...
        % 'MSH_posterior_1min_mean', 'SW_posterior_1min_mean', 'MSP_posterior_1min_mean', ...
        % 'MSH_posterior_2min_mean', 'SW_posterior_2min_mean', 'MSP_posterior_2min_mean', ...
        % 'MSH_posterior_4min_mean', 'SW_posterior_4min_mean', 'MSP_posterior_4min_mean', ...
        % 'MSH_posterior_8min_mean', 'SW_posterior_8min_mean', 'MSP_posterior_8min_mean', ...
        % ...
        % 'MSH_mahaldist_15sec_mean', 'SW_mahaldist_15sec_mean', 'MSP_mahaldist_15sec_mean', ...
        % 'MSH_mahaldist_30sec_mean', 'SW_mahaldist_30sec_mean', 'MSP_mahaldist_30sec_mean', ...
        % 'MSH_mahaldist_1min_mean', 'SW_mahaldist_1min_mean', 'MSP_mahaldist_1min_mean', ...
        % 'MSH_mahaldist_2min_mean', 'SW_mahaldist_2min_mean', 'MSP_mahaldist_2min_mean', ...
        % 'MSH_mahaldist_4min_mean', 'SW_mahaldist_4min_mean', 'MSP_mahaldist_4min_mean', ...
        % 'MSH_mahaldist_8min_mean', 'SW_mahaldist_8min_mean', 'MSP_mahaldist_8min_mean', ...
        % ...
        % 'logpdf_15sec_mean', ...
        % 'logpdf_30sec_mean', ...
        % 'logpdf_1min_mean', ...
        % 'logpdf_2min_mean', ...
        % 'logpdf_4min_mean', ...
        % 'logpdf_8min_mean'
    % };

    n_cols = length(header);
    column_id = split(num2str(1:n_cols));
    parenthesis = join([repelem({'('}, n_cols)', column_id, repelem({')'}, n_cols)'], '');
    header = join(join([parenthesis, header'], ''), ', ');
    header = header{:};
    header = [header, '\n'];
end