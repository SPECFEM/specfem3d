%
% function [m_str, m_ts, m_xs, m_ys] = wave2d_splitm(m,m_inds)
% Carl Tape, 25-Jan-2010
%
% This function splits a wave2d.f90 model vector into constituent parts.
%
% calls xxx
% called by xxx
%

function [m_str, m_ts, m_xs, m_ys] = wave2d_splitm(m,m_inds)

if length(m) ~= m_inds(4,2)
    whos m
    m_inds
    error('incompatible indexing');
end

m_str = m(m_inds(1,1):m_inds(1,2));
m_ts  = m(m_inds(2,1):m_inds(2,2));
m_xs  = m(m_inds(3,1):m_inds(3,2));
m_ys  = m(m_inds(4,1):m_inds(4,2));

% ensure that the output are columns
m_str = m_str(:);
m_ts  = m_ts(:);
m_xs  = m_xs(:);
m_ys  = m_ys(:);

%=========================================================