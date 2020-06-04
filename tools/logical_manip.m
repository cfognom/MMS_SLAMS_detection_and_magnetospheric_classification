function logic = logical_manip(logic, method)
%LOGICAL_MANIP Summary of this function goes here
%   Detailed explanation goes here
    % logical(logic)
    if nargin == 1
        method = 'firstOfSimilar';
    end
    dif = diff(logic);
    switch method
        case 'firstOfSimilar'
            logic = [true, dif ~= 0];
        case 'firstOfOne'
            logic = [logic(1) == 1, dif == 1];
        case 'firstOfZero'
            logic = [logic(1) == 0, dif == -1];
        case 'lastOfSimilar'
            logic = [dif ~= 0, true];
        case 'lastOfOne'
            logic = [dif == -1, logic(end) == 1];
        case 'lastOfZero'
            logic = [dif == 1, logic(end) == 0];
    end
end

