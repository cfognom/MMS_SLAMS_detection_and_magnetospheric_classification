function BS = GSE2BS(GSE)
%GSE2BS Summary of this function goes here
%   Detailed explanation goes here
    xy = irf_abs(GSE(:,2:3), 1);
    xBS = GSE(:,1) - bowshock_pos(xy);
    BS = [xBS, xy];
end

