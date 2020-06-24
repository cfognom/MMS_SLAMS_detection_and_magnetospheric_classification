function testfunc()
%TESTFUNC Summary of this function goes here
%   Detailed explanation goes here

    fileID = fopen('test_mat.txt', 'r');
    % fgetl(fileID)
    tmp = fscanf(fileID, '%d,%d = [%f%f%f%f]', [2,10])
    fclose('all');
    % tmp(:,1)
    % tmp(:,2)
end

