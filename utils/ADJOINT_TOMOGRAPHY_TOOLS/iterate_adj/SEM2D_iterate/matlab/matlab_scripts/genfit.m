%
% function [m1, e1] = genfit(theory, m0, jogs, data, param)
% Carl Tape, 23-Sep-2002
% printed 23-Sep-2002
%
% Program by John Woodhouse, Dept. Earth Sciences, University of Oxford.
% Non-linear least squares algorithm to fit data
% to a model using the theory function 'theory'.
%
% original Oxford location =  /home/owl_data2/john/matlab/genfit.m
%
% theory  = a function that generates a column
%           of predicted data given a model m and 
%           having the call: d = theory(m, param).
% m0      = starting model
% jogs    = a column of perturbations to the
%           model that genfit will use to estimate the
%           partial derivatives of theory with respect to
%           the model parameters. A zero value in jogs
%           indicates that the corresponding model parameter
%           is to be held fixed.
% data    = vector of data
% param   = an array of parameters associated with the data,
%           having the same number of rows as data
%          
% m1      = the updated model.
% e1      = an estimate of the uncertainties (standard errors)
%           of the updated model.
%
% calls theory.m (e.g., theoryHyp.m)
% called by testgenfitA.m, testgenfitB.m, testgenfitC.m
%

function [m1, e1] = genfit(theory, m0, jogs, data, param)

numd = length(data);   % number of data points
numm = length(m0);     % number of model parameters

resid = data - feval(theory, m0, param);
r0 = resid' * resid;
rms_residual = sqrt(r0/numd);

%disp('  ');
%disp('  data    theory    resid');
%disp([data feval(theory, m0, param) resid]);
    
% vary each model parameter in turn, to calculate partial derivatives
a  = zeros(numd, numm);

% predicted data from model m0
f0 = feval(theory, m0, param);

% Calculate partial derivatives of 'theory' with respect
% to each model parameter, storing the results in matrix a.
for i=1:numm
    xx = m0;
    dx = jogs(i);
    if dx ~= 0,                     % check for non-vanishing dx
        xx(i) = xx(i) + dx;         % new model
        a(:, i) = (feval(theory, xx, param) - f0) ./ dx;
    else
        a(:, i) = zeros(numd, 1);
    end
end

% Use the 'pseudo-inverse' pinv(a'*a), rather than the
% true inverse inv(a'*a), because a'a may have zero rows and columns.
c1 = pinv(a' * a);
m1 = m0 + c1*(a' * resid);          % updated model

c1 = c1 * rms_residual^2;
e1 = zeros(numm, 1);
for i = 1:numm,
    e1(i) = sqrt(c1(i, i));         % formal standard errors in model parameters
end

%==========================================================
