function y = classifier(X)
%CLASSIFIER Summary of this function goes here
%   Detailed explanation goes here
    [n, c] = size(X);
    y = zeros(n,1);
    for i = 1:n
        if X(i,3) < (-0.01*X(i,2) + 5)
            y(i) = 1;
        else
            % if X(i,3) < (0.4*X(i,2) - 100)
            if X(i,3) < (0.25*X(i,2) - 60) && X(i,3) < (-0.035*X(i,2) + 28) % && X(i,1) < (0.1*X(i,2) - 22) && X(i,1) < (-0.01*X(i,2) + 13)
                y(i) = 2;
            else
                y(i) = 3;
            end
        end
    end
end

