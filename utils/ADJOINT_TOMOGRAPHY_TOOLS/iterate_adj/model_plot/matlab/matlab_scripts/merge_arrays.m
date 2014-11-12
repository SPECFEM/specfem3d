%
% function out_array = merge_arrays(array1,array2,stconnect)
% CARL TAPE, 08-Aug-2007
% printed xxx
%
% This file inputs two arrays of strings and merges them with a connecting
% character.
%
% See also split_array.m, which is the opposite.
%
% calls xxx
% called by xxx
%

function out_array = merge_arrays(array1,array2,stconnect)

num = length(array1);
if length(array1) ~= length(array2)
    error(' array1 and array2 must have the same length');
end

% allocate array; initialize array
% THIS SPEEDS THINGS UP BY A FACTOR OF ABOUT 100.
out_array = cellstr(repmat(' ',num,1));

for ii = 1:num
    %if mod(ii,500)==0, disp([num2str(ii) ' out of ' num2str(num)]); end
    out_array{ii} = [array1{ii} stconnect array2{ii}];
end
out_array = out_array(:);

%=======================================================================