function suffix=suf(string,delim,iter)
% suffix=SUF(string,delim,iter)
% 
% Finds the suffix or extension in a string, say a filename.
%
% INPUT:
%
% string        The string to be parsed
% delim         The string that separates prefix from suffix
% iter          0 Don't iterate
%               1 Iterate on the suffix to get the last bit
% 
% See also PREF
%
% Last modified by fjsimons-at-alum.mit.edu, 05/06/2008

defval('delim','.');
defval('iter',1)

ldel=length(delim);
[prefix,suffix]=strtok(string,delim);

% This added 17.11.2004
if iscell(suffix)
  suffix=cell2mat(suffix);
end

suffix=suffix(ldel+1:end);

if iter==1
  % This suffix might not be the last one
  while findstr(suffix,delim)
    suffix=suf(suffix,delim);
  end
end




