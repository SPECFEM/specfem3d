%
% written by Carl Tape, University of Fairbanks, Alaska, USA
%
function [Mout,T] = convert_MT(i1,i2,M)
%CONVERT_MT convert moment tensor matrices among different bases
%
% This program converts between different moment tensor conventions.
% All conventions are associated with a local coordinate system.
%
% M = [M11 M22 M33 M12 M13 M23]
%
% INPUT
%   i1      index of input moment tensor basis
%   i2      index of output moment tensor basis
%   M       6 x n set of moment tensors, M = [M11 M22 M33 M12 M13 M23]
%
% OUTPUT
%   Mout    6 x n set of moment tensors in basis of i2
%   T       transformation matrix to change basis of M from i1 to i2: Mout = T*M*T'
%
% Convention 1: up-south-east (GCMT) (www.globalcmt.org), also CMTSOLUTION format
%   1: up (r), 2: south (theta), 3: east (phi)
%
% Convention 2: Aki and Richards (1980, p. 114-115, 118)
%   also Jost and Herrman (1989, Fig. 1)
%   1: north, 2: east, 3: down
%
% Convention 3: Stein and Wysession (2003, p. 218)
%   also TapeTape2012a "A geometric setting for moment tensors" (p. 478)
%   also several Kanamori codes
%   1: north, 2: west, 3: up
% 
% Convention 4: 
%   1: east, 2: north, 3: up
% 
% Convention 5: TapeTape2013 "The classical model for moment tensors" (p. 1704)
%   1: south, 2: east, 3: up
%
% See the vector version of this function: convertv.m 
%
% Carl Tape, 11/2010
%

NTYPE = 5;  % number of right-handed bases to consider

% permutation matrices
Tall = cell(NTYPE,NTYPE);
Tall{1,1} = eye(3);
Tall{2,2} = eye(3);
Tall{3,3} = eye(3);
Tall{4,4} = eye(3);
Tall{5,5} = eye(3);
% from i1 to i2
Tall{1,2} = [0 -1 0 ; 0 0 1 ; -1 0 0];
Tall{1,3} = [0 -1 0 ; 0 0 -1 ; 1 0 0];
Tall{1,4} = [0 0 1 ; 0 -1 0 ; 1 0 0];
Tall{1,5} = [0 1 0 ; 0 0 1 ; 1 0 0];
Tall{2,3} = [1 0 0 ; 0 -1 0 ; 0 0 -1];
Tall{2,4} = [0 1 0 ; 1 0 0 ; 0 0 -1];
Tall{2,5} = [-1 0 0 ; 0 1 0 ; 0 0 -1];
Tall{3,4} = [0 -1 0 ; 1 0 0 ; 0 0 1];
Tall{3,5} = [-1 0 0 ; 0 -1 0 ; 0 0 1];
Tall{4,5} = [0 -1 0 ; 1 0 0 ; 0 0 1];
% from i2 to i1
Tall{2,1} = Tall{1,2}';
Tall{3,1} = Tall{1,3}';
Tall{3,2} = Tall{2,3}';
Tall{4,1} = Tall{1,4}';
Tall{4,2} = Tall{2,4}';
Tall{4,3} = Tall{3,4}';
Tall{5,1} = Tall{1,5}';
Tall{5,2} = Tall{2,5}';
Tall{5,3} = Tall{3,5}';
Tall{5,4} = Tall{4,5}';
% transformation matrix
T = Tall{i1,i2};

stlabs = {'up-south-east (GCMT)',...
    'north-east-down (AkiRichards)',...
    'north-west-up',...
    'east-north-up',...
    'south-east-up'};
disp(sprintf('convert_MT.m: %s to %s',stlabs{i1},stlabs{i2}));

if nargin==2
   disp('returning transformation matrix only, as requested');
   Mout = T;
   return
end

% make sure M is 6 x n
[M,n] = Mdim(M);

if i1==i2
    %error('i1 must differ from i2');
    disp('warning: i1 = i2, so no change');
    Mout = M;
    return
end

Mout = [];  % initialize

if i1==1
    if i2==2        % up-south-east (GCMT) to north-east-down (AkiRichards) (AR, 1980, p. 118)
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(3,:);
        Mout(3,:) = M(1,:);
        Mout(4,:) = -M(6,:);
        Mout(5,:) = M(4,:);
        Mout(6,:) = -M(5,:);
    elseif i2==3    % up-south-east (GCMT) to north-west-up (/opt/seismo-util/bin/faultpar2cmtsol.pl)
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(3,:);
        Mout(3,:) = M(1,:);
        Mout(4,:) = M(6,:);
        Mout(5,:) = -M(4,:);
        Mout(6,:) = -M(5,:);
    elseif i2==4    % up-south-east (GCMT) to east-north-up
        Mout(1,:) = M(3,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(1,:);
        Mout(4,:) = -M(6,:);
        Mout(5,:) = M(5,:);
        Mout(6,:) = -M(4,:);
    elseif i2==5    % up-south-east (GCMT) to south-east-up
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(3,:);
        Mout(3,:) = M(1,:);
        Mout(4,:) = M(6,:);
        Mout(5,:) = M(4,:);
        Mout(6,:) = M(5,:);  
    end
    
elseif i1==2
    if i2==1        % north-east-down (AkiRichards) to up-south-east (GCMT) (AR, 1980, p. 118)
        Mout(1,:) = M(3,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(2,:);
        Mout(4,:) = M(5,:);
        Mout(5,:) = -M(6,:);
        Mout(6,:) = -M(4,:);
    elseif i2==3    % north-east-down (AkiRichards) to north-west-up
        Mout(1,:) = M(1,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = -M(5,:);
        Mout(6,:) = M(6,:);   
    elseif i2==4    % north-east-down (AkiRichards) to east-north-up
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = M(4,:);
        Mout(5,:) = -M(6,:);
        Mout(6,:) = -M(5,:);
    elseif i2==5    % north-east-down (AkiRichards) to south-east-up
        Mout(1,:) = M(1,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = M(5,:);
        Mout(6,:) = -M(6,:);   
    end
    
elseif i1==3
    if i2==1        % north-west-up to up-south-east (GCMT)
        Mout(1,:) = M(3,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(2,:);
        Mout(4,:) = -M(5,:);
        Mout(5,:) = -M(6,:);
        Mout(6,:) = M(4,:);
    elseif i2==2    % north-west-up to north-east-down (AkiRichards)
        Mout(1,:) = M(1,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = -M(5,:);
        Mout(6,:) = M(6,:); 
    elseif i2==4    % north-west-up to east-north-up
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = -M(6,:);
        Mout(6,:) = M(5,:); 
    elseif i2==5    % north-west-up to south-east-up
        Mout(1,:) = M(1,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = M(4,:);
        Mout(5,:) = -M(5,:);
        Mout(6,:) = -M(6,:); 
    end
    
elseif i1==4
    if i2==1        % east-north-up to up-south-east (GCMT)
        Mout(1,:) = M(3,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(1,:);
        Mout(4,:) = -M(6,:);
        Mout(5,:) = M(5,:);
        Mout(6,:) = -M(4,:);
    elseif i2==2    % east-north-up to north-east-down (AkiRichards)
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = M(4,:);
        Mout(5,:) = -M(6,:);
        Mout(6,:) = -M(5,:);
    elseif i2==3    % east-north-up to north-west-up
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = M(6,:);
        Mout(6,:) = -M(5,:); 
    elseif i2==5    % east-north-up to south-east-up
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = -M(6,:);
        Mout(6,:) = M(5,:); 
    end
    
elseif i1==5        % south-east-up to up-south-east (GCMT)
    if i2==1
        Mout(1,:) = M(3,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(2,:);
        Mout(4,:) = M(5,:);
        Mout(5,:) = M(6,:);
        Mout(6,:) = M(4,:);
    elseif i2==2    % south-east-up to north-east-down (AkiRichards)
        Mout(1,:) = M(1,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = M(5,:);
        Mout(6,:) = -M(6,:);
    elseif i2==3    % south-east-up to north-west-up
        Mout(1,:) = M(1,:);
        Mout(2,:) = M(2,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = M(4,:);
        Mout(5,:) = -M(5,:);
        Mout(6,:) = -M(6,:);
    elseif i2==4    % south-east-up to east-north-up
        Mout(1,:) = M(2,:);
        Mout(2,:) = M(1,:);
        Mout(3,:) = M(3,:);
        Mout(4,:) = -M(4,:);
        Mout(5,:) = M(6,:);
        Mout(6,:) = -M(5,:); 
    end
    
end

%==========================================================================
% EXAMPLES

if 0==1
    % transformation matrix only
    i1 = 1; i2 = 2;
    T = convert_MT(i1,i2)
    
    % simple example
    i1 = 1; i2 = 2;
    A = [1:6]'
    M1 = convert_MT(i1,i2,A)
    M2 = convert_MT(i2,i1,M1)
    
    % checking the transformation matrix
    i1 = 1; i2 = 5;
    M1 = rand(6,1);
    [M2,T] = convert_MT(i1,i2,M1);
    Mvec2Mmat(M1,1)                 % up-south-east
    Mcheck = T*Mvec2Mmat(M1,1)*T'   % south-east-up
    Mvec2Mmat(M2,1)                 % (check)
    % example vector v, example matrix M
    v1 = rand(3,1)                  % up-south-east
    v2 = T*v1                       % south-east-up
    X = randi(10,3); X1=X'*X,  X2=T*X1*T'
    
    % check all possible transofrmations
    M1 = rand(6,1);
    NTYPE = 5;
    for i1=1:NTYPE
        for i2=2:NTYPE
            if i2==i1, continue; end
            [M2,T] = convert_MT(i1,i2,M1);
            if and(~isempty(T),~isempty(M2))
                Mcheck = T*Mvec2Mmat(M1,1)*T';
                ncheck = norm(Mvec2Mmat(M2,1)-Mcheck);
                % display info if the numerical check fails
                if ncheck > 1e-6
                   disp(sprintf('from %i to %i',i1,i2));
                   Mvec2Mmat(M2,1),Mcheck
                   v1 = rand(3,1) 
                   v2 = T*v1
                   error('check') 
                end
            end
        end
    end
end

%==========================================================================

