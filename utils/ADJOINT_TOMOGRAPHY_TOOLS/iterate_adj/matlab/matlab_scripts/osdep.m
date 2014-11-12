function of=osdep
% of=OSDEP
%
% Returns the value of the read option to be used
% on the LOCAL operating system for files created
% on the LOCAL operating system.
%
% Last modified by fjsimons-at-alum.mit.edu, 23.11.2004


if strcmp(getenv('OSTYPE'),'linux')
  of= 'l';  
end
if strcmp(getenv('OSTYPE'),'solaris')
  of= 'b';  
end
