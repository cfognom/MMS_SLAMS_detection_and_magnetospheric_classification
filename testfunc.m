function testfunc()
%TESTFUNC Summary of this function goes here
%   Detailed explanation goes here

    fileID = fopen('test_cell.csv', 'r');
    fgetl(fileID)
    tmp = textscan(fileID, '%d %s %f %f %f', 'Delimiter', ',')
    fclose('all');
    tmp{:,1}
    tmp{:,2}
    tmp{:,3}
    tmp{:,4}
    tmp{:,5}
end

