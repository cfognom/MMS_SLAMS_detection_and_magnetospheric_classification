function logic = logical_manip(logic, method)
%LOGICAL_MANIP Summary of this function goes here
%   Detailed explanation goes here
    % logical(logic)
    if nargin == 1
        method = 'firstOfSimilar';
    end
    [row, col] = size(logic);
    if row == 1
        logic = execute(logic);
    elseif col == 1
        logic = execute(logic');
        logic = logic';
    else
        error('logical_manip can only process vectors. Either n_row or n_col must be 1.');
    end

    function logic = execute(logic) 
        switch method
            case 'firstOfSimilar'
                logic = [true, xor(logic(2:end), logic(1:end - 1))];
            case 'firstOfOne'
                logic = [logic(1) == 1, logic(2:end) & ~logic(1:end - 1)];
            case 'firstOfZero'
                logic = [logic(1) == 0, ~logic(2:end) & logic(1:end - 1)];
            case 'lastOfSimilar'
                logic = [xor(logic(2:end), logic(1:end - 1)), true];
            case 'lastOfOne'
                logic = [~logic(2:end) & logic(1:end - 1), logic(end) == 1];
            case 'lastOfZero'
                logic = [logic(2:end) & ~logic(1:end - 1), logic(end) == 0];
        end
    end
end

