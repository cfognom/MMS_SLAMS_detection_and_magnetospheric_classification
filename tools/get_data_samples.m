function X = get_data_samples(filePrefix, variable, t_samples, time_window, varargin)
%GET_DATA_SAMPLES Summary of this function goes here
%   Samples 'variable' from 'filePrefix' at 't_samples' which is an EpochTT array of time instants.
%   'time_threshold' is how much extra time in seconds should be loaded before and after each point.
%   Name-value pairs:
%       'modFunc', function handle to a function that will be applied to the data of the TSeries before sampling.


    p = inputParser;
    addParameter(p, 'modFunc', @linear_interpolation);
    parse(p, varargin{:});
    r = p.Results;

    % Init
    t_samples = sort(t_samples);
    n_samples = length(t_samples);
    % n_instruments = length(parameters);
    % n_variables = cellfun(@(x) length(x{3}), parameters);
    % n_variables_tot = sum(n_variables);

    % % Get dt for each instrument, to be used later to make a small interval for loading
    % dt = zeros(n_instruments,1);
    % for i = 1:n_instruments
    %     p = parameters{i};
    %     ts = load_tints_MMS(p{1}, p{2}(1), [t_samples(1), t_samples(1)]);
    %     t_res = ts.userData.GlobalAttributes.Time_resolution{1};
    %     longest_delay = max(sscanf(t_res(~isletter(t_res)), '%f'));
    %     dt(i) = 1.2*longest_delay; % Safety factor of 1.2
    % end

    dt = time_window/2;
    done = zeros(n_samples, 1);
    % loaded = zeros(n_samples, 1);
    X = zeros(n_samples, 1);
    while any(~done)
        t_not_done = t_samples(~done);
        t = t_not_done(1);
        tint_load = [t + -dt, t + dt];
        ts = load_tints_MMS(filePrefix, variable, tint_load, 'trim', false);
        loaded = t_samples < (ts.time(end) + -dt);
        batch = (~done & loaded);
        if ~any(batch)
            loaded_partial = t_samples < ts.time(end);
            batch = (~done & loaded_partial);
        end
        first = find(batch, 1, 'first');
        last = find(batch, 1, 'last');
        t_batch = t_samples(batch);
        data = r.modFunc(ts, t_batch, dt);
        [~, d] = size(data);
        X(first:last, 1:d) = data;
        done = done | batch;
        fprintf('Getting data samples for variable ''%s'' %u/%u\n', variable, nnz(done), n_samples);
    end

    function out = linear_interpolation(ts, t, ~)
        out = interp1(double(ts.time.epoch), ts.data, double(t.epoch));
    end

    % % Load small interval for each point, interpolate and place into X
    % X = zeros(n_trainPoints, n_variables_tot);
    % bad_points = zeros(1, n_trainPoints);
    % for tp = 1:n_trainPoints
    %     fprintf('Getting train data %u/%u\n', tp, n_trainPoints);
    %     c = 1;
    %     t = double(t_samples(tp).epoch);
    %     for i = 1:n_instruments
    %         p = parameters{i};
    %         tint = [t_samples(tp) + -dt(i), t_samples(tp) + dt(i)];
    %         try
    %             TS = MMS_load(p{1}, p{2}, tint);
    %             % TS = TSeries.empty;
    %             for v = 1:n_variables(i)
    %                 if iscell(TS)
    %                     ts = TS{v};
    %                 else
    %                     ts = TS;
    %                 end
    %                 if length(p) > 2 && p{3}(v) == "abs"
    %                     ts = irf_abs(ts);
    %                 end
    %                 [~, d] = size(ts.data);
    %                 X(tp,c:c + d - 1) = interp1(double(ts.time.epoch), ts.data, t, 'linear', 'extrap');
    %                 c = c + d;
    %             end 
    %         catch
    %             bad_points(tp) = 1;
    %         end
    %     end
    % end

    % % Remove bad points, i.e. points that are missing some instrument data.
    % X = X(bad_points == 0,:);
end

