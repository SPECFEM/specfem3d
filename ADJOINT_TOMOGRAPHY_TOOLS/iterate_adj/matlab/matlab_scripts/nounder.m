function sin=nounder(sin)
% sin=NOUNDER(sin)
%
% Changes underscores to dashes in a string
%
% Last modified by fjsimons-at-alum.mit.edu, Feb 06th, 2004

sin(find(abs(sin)==95))='-';

