function x = bowshock_pos(yz)
%BOWSHOCK_POS Summary of this function goes here
%   Detailed explanation goes here
    K = 25;
    epsilon = 0.8;
    x = 2*(epsilon*K - sqrt(K^2 + 0.25*(4*epsilon^2 - 4)*yz.^2))/(2*epsilon^2 - 2);
end

