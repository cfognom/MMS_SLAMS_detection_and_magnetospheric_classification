function ts = load_tints_MMS(filePrefix, variable, tints, varargin)
    %MMS_LOAD loads MMS data
    %   IN:
    %       filePrefix: String of the name of the file of interest without date.
    %       variable: String of the variable of interest.
    %       tints: EpochTT time intervals.
    %   Out:
    %       ts: TSeries

    p = inputParser;
    addParameter(p, 'trim', true);
    parse(p, varargin{:});
    r = p.Results;

    n_tints = length(tints)/2;

    for i = 1:n_tints
        tint = select_tint(tints, i);
        if r.trim
            tmp_ts = mms.db_get_ts(filePrefix, variable, tint);
            if iscell(tmp_ts)
                tmp_tmp_ts = tmp_ts{1};
                for j = 2:length(tmp_ts)
                    tmp_tmp_ts = tmp_tmp_ts.combine(tmp_ts{j});
                end
                tmp_ts = tmp_tmp_ts;
            end
        else
            tmp_ts = mms.variable2ts(mms.db_get_variable(filePrefix, variable, tint));
        end
        if i == 1
            ts = tmp_ts;
        else
            ts = ts.combine(tmp_ts);
        end
    end

    % function ts = getTS(filePrefix, variable, tint)
    %     ts = mms.db_get_ts(char(filePrefix), char(variable), tint);
    %     if class(ts) == "double"
    %         error("Loading variable '%s' from file type '%s' during the time interval: '%s' to '%s' returned an empty time series.", sprintf(varName, variable), char(filePrefix), tint(1).utc, tint(2).utc);
    %     end
    % end
end

