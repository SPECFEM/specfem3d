%
% function xmin = quad_min_4(P,x0)
% Carl Tape, 11-Jan-2006
% printed xxx
%
% This function inputs two (x,y) points and one slope, and returns a
% quadratic fitting polynomial, along with the (analytical) minimum value.
%
% calls quad_shift.m
% called by model_optimize.m
%

function [xmin,P1] = quad_min_4(x1,x2,y1,y2,g1,opts,stlabs)

ifig = opts(1);     % =1 to plot figure
isub = opts(2);     % =1 to plot as a subfigure; =0 to plot as a full figure

if length([x1 x2 y1 y2 g1]) ~= 5, error('check input dimensions'); end

a = ((y2 - y1) - g1*(x2 - x1)) / (x2^2 - x1^2);
b = g1;
c = y1 - a*x1^2 - b*x1;
%c = y2 - a*x2^2 - b*x2;

% ax^2 + bx + c
P1 = [a b c]';

% a(x-b)^2 + c
[P2,qvert,stit] = quad_shift(P1,1);
xmin = qvert(1);

if ifig==1
    if isub==1
        specs = [1 6 12 10]; fac = 0.1;
    else
        figure;
        specs = [2 14 18 12]; fac = 0.05;
    end
    
    % step COULD be negative (source inversion)
    temp = sort([x1 x2]);
    x1plot = temp(1);
    x2plot = temp(2);
    if x1 ~= x1plot
        iflip = 1;
        stlabs = [stlabs(1:2) fliplr(stlabs(3:5))]
    end
    
    axpoly = axes_expand([x1plot x2plot 0 max([y1 y2])],1.2,1);
    axpoly(3) = 0;
    dy = axpoly(4) - axpoly(3);
    dx = axpoly(2) - axpoly(1);
    ylab = axpoly(3) - fac*dy;
    
    ymin = polyval(P1,xmin);  % quadratic function evaluation

    % base level for test-model parabola
    K = 0.5;
    %K = 0.0;   % no model norm term or data errors
    
    % initial guess is based on a quadratic fit
    %aquad = g1^2/(4*y1);
    aquad = -g1^2/(4*(K - y1));
    Pquad = [aquad g1 y1]';
    
    % x-points for smooth curves
    n = 100; xpts = linspace(axpoly(1),axpoly(2),n);
    
    % curves through (x1,y1)
    g1_line = polyval([g1 y1-g1*x1],xpts);
    g1_test = polyval(Pquad,xpts);
    g1_quad = polyval(P1,xpts);

    hold on;
    
    % plot curves
    plot(xpts,g1_quad,'b','linewidth',specs(1));
    plot(xpts,g1_test,'b--');
    plot(xpts,g1_line,'r--');
    
    % plot black lines
    plot([x1 x1 axpoly(1)],[0 y1 y1],'k');
    plot([x2 x2 axpoly(1)],[0 y2 y2],'k');
    plot([xmin xmin axpoly(1)],[0 ymin ymin],'k--');
    
    % plot markers
    plot([x1 x2],[y1 y2],'ko','markersize',specs(2),'MarkerFaceColor','b');
    plot(xmin,ymin,'bo',x2,0,'bo',xmin,0,'bo','markersize',specs(2),'MarkerFaceColor','w');
    %plot(0.5*x2,0,'go','markersize',specs(2),'MarkerFaceColor','w');  % linear projection

    axis(axpoly);
    xlabel(stlabs{1},'fontsize',specs(3));
    ylabel(stlabs{2},'fontsize',specs(3));
    if isub==0
        title({stit{1},stit{2}},'fontsize',specs(3))
        grid on;
    else
        set(gca,'xtick',[x1plot x2plot],'xticklabel',{[],[]}); 
        %set(gca,'xtick',[x1plot xmin x2plot]),'xticklabel',{'0',[],[]};
    end
    text(x1plot,ylab,stlabs{3},'fontsize',specs(4));
    text(xmin,ylab,stlabs{4},'fontsize',specs(4));
    text(x2plot,ylab,stlabs{5},'fontsize',specs(4));
     orient tall, wysiwyg
end

if 0==1
   x1 = 0
   x2 = randomint(1,5,1)
   y1 = randomint(1,5,1)
   y2 = randomint(1,5,1)
   g1 = randomint(-5,-1,1)

   opts = [1 0];
   stk  = num2str(0);
   stlabs = {'\lambda',['\chi^{' stk '} ( \lambda )'],'0',['\lambda_{' stk '}'],['\lambda_{' stk 't}']};
   
   [xmin,P1] = quad_min_4(x1,x2,y1,y2,g1,opts,stlabs);
   set(gca,'xtick',[-10:10],'ytick',[-10:10]);
   axis equal;
end

%=========================================================