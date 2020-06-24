function tints = UTC2tints(UTC)
    UTC = UTC';
    UTC = reshape(UTC, [], 1);
    UTC = char(UTC);
    tints = EpochTT(UTC);
end

