function x = bowshock_pos(yz)
%BOWSHOCK_POS Summary of this function goes here
%   Detailed explanation goes here
    K = 25;
    epsilon = 0.8;
    x = (epsilon*K - sqrt(K^2 + (epsilon^2 - 1)*yz.^2))/(epsilon^2 - 1);
end

