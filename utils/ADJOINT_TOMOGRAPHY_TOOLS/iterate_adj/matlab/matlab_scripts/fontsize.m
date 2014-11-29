%
% function function fontsize(fs,h)
% Juliette Artru, 15-Sept-2004
% printed xxxx
%
% This changes the font on an entire figure.
%
% fs: font size
% h: handle graphics (optional)
%

function fontsize(fs,h);

if(nargin<2);h=gcf;end
hc=get(gcf,'child');
hall=hc;
for k=1:length(hc);
    if(strcmp(get(hc(k),'Type'),'axes'))
        hall=[hall; get(hc(k),'XLabel') ; get(hc(k),'YLabel') ; get(hc(k),'Title')];
    end
end
set(hall,'fontsize',fs);
  
 
