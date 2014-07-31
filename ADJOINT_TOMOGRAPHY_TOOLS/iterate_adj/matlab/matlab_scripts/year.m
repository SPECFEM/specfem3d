function y = year(d) 
%YEAR Year of date. 
%   Y = YEAR(D) returns the year of a serial date number or a date string, D. 
% 
%   For example, y = year(728647) or y = year('19-Dec-1994') returns y = 1994. 
%  
%   See also DATEVEC, DAY, MONTH. 
 
%   Copyright 1995-2006 The MathWorks, Inc.  
%   $Revision: 1.6.2.2 $   $Date: 2006/04/20 17:47:27 $ 
 
if nargin < 1 
  error('finance:year:missingInputs','Please enter D.') 
end 
if ischar(d) 
  d = datenum(d); 
end 
 
c = datevec(d(:));      % Generate date vectors from dates 
y = c(:,1);             % Extract years  
if ~ischar(d) 
  y = reshape(y,size(d)); 
end 

