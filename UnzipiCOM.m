function [outputFile] = UnizipiCOM(fn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% for 7zip 15.14 x64 
file7z = fn;
fn_chk = fn(1,1:end-3)
if isfile(fn_chk)
    return
end
% rest of code
[dir7z,~,~] = fileparts( file7z );
dir7z = ['"' dir7z '"'];
file7z = ['"' file7z '"'];
[status,cmdout] = system( ['"C:\Program Files\7-Zip\7z.exe" x -o',dir7z,' ',file7z] )

end

