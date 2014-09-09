%
% function mvec = wave2d_m2mvec(m,m_inds,beta0)
% Carl Tape, 25-Jan-2010
%
% Convert a model vector from 'inversion parameters' to 'physical
% parameters'.  In our case, only the structural parameters differ.
%
% calls xxx
% called by xxx
%

function mvec = wave2d_m2mvec(m,m_inds,beta0)

nmod = length(m);
nmod_str = m_inds(1,2);

mvec = zeros(nmod,1);
mvec(1:nmod_str) = beta0 * exp( m(1:nmod_str) );
mvec(nmod_str+1 : nmod) = m(nmod_str+1 : nmod);

%=========================================================
