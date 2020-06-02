function epochtt = datetimeToEpochTT(datetime)
%DATETIMETOEPOCHTT Summary of this function goes here
%   Detailed explanation goes here
    % [row, col] = size(datetime);
    % % epochtt(row, col) = EpochTT();
    % for i = row:-1:1
    %     for j = col:-1:1
    %         epochtt(i,j) = irf_time(datevec(datetime(i,j)), 'vector>epochtt');     
    %     end
    % end
    epochtt = EpochTT(datestr(datetime, 'yyyy-mm-ddTHH:MM:SS.FFFZ'));
end

