function defval(name,value)
% DEFVAL(name,value)
%
% Checks if a variable 'name' exists in the caller's workspace,
% if nonexistent or if empty 'name' is set to value 'value' in caller's workspace
%
% Won't work for an unassigned structure variable
% 
% Last modified by fjsimons-at-alum.mit.edu, 09/12/2007

if ~ischar(name),
  error('The first argument of defval has to be a string (variable name)');
end

si=1;
if evalin('caller',[ 'exist(''' name ''')']);,
  si=evalin('caller',[ 'isempty(' name ')']);
end
if si,
  assignin('caller',name,value);
  na=dbstack;
  %disp(['Default value used in ' na(2).name ': ' name '=' num2str(value(1,:))])
end
  
