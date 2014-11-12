%
% function [Pout,qvert,stit] = quad_shift(Pin,opts)
% Carl Tape, 11-Jan-2006
% printed xxx
%
% This function adjusts parabola (i.e., quadratic polynomial) coefficients 
% between two representations, and also returns the vertex.
%
% See also cubic_shift.m
%
% calls xxx
% called by test_poly.m
%

function [Pout,qvert,stit] = quad_shift(Pin,opts)

iopt = opts(1);
stx = '%.3f';

P = Pin(:);
a = Pin(1);
b = Pin(2);
c = Pin(3);

% adjust the polynomial coefficients
if iopt == 1
    % ax^2 + bx + c --> a(x-b)^2 + c
    Pout(1) = a;
    Pout(2) = -b/(2*a);
    Pout(3) = c - b^2/(4*a);
    
    qvert = [Pout(2) Pout(3)];
    
    stit1 = ['y = ' num2str(sprintf(stx, a)) ' x^2  +  ' ...
                    num2str(sprintf(stx, b)) ' x  +  ' ...
                    num2str(sprintf(stx, c)) ];  
    stit2 = ['y = ' num2str(sprintf(stx, Pout(1))) ' (x - ' ...
                    num2str(sprintf(stx, Pout(2))) ')^2  +  ' ...
                    num2str(sprintf(stx, Pout(3))) ];
    
else
    % a(x-b)^2 + c --> ax^2 + bx + c
    if a==0
        Pout(1) = 0;
        Pout(2) = 0;
        Pout(3) = c;
        qvert = [NaN NaN];
    else
        Pout(1) = a;
        Pout(2) = -2*a*b;
        Pout(3) = a*b^2 + c;
        qvert = [b c];
    end
    
    stit1 = ['y = ' num2str(sprintf(stx, Pout(1))) ' x^2  +  ' ...
                    num2str(sprintf(stx, Pout(2))) ' x  +  ' ...
                    num2str(sprintf(stx, Pout(3))) ];
    stit2 = ['y = ' num2str(sprintf(stx, a)) ' (x - ' ...
                    num2str(sprintf(stx, b)) ')^2  +  ' ...
                    num2str(sprintf(stx, c)) ];
end

Pout = Pout(:);
stit = {stit1,stit2};

%=========================================================
