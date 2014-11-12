%
% function norm_parts = wave2d_norm_sq(m,mprior,C,m_inds,weights,imnorm)
% Carl Tape, 25-Jan-2010
%
% Compute the norms of different parts of a model vector or gradient.
%
% calls xxx
% called by xxx
%

function norm_parts = wave2d_norm_sq(m,mprior,C,m_inds,weights,imnorm)

m = m(:);
[a,b] = size(C);
if a==1, error('check dimension of C'); end
if b==1, icdiag = 1; else icdiag = 0; end

npar = length(m_inds);

%stfmt = repmat('%16.4e',1,npar);

Cmat = zeros(a,b);
mvec = zeros(length(m),1);
npw = zeros(npar,1);
if imnorm==0
    Cmat = C;
    mvec = m;
    npw = 1 ./ weights;
    
elseif imnorm==1
    %if icdiag==1, Cmat = 1./C; else Cmat = inv(C); end
    Cmat = diag(1./diag(C));  % APPROXIMATION
    mvec = m - mprior;
    npw = weights;
end

disp('norm of parts:');
norm_parts = zeros(npar,1);
for ii=1:npar
    inds = m_inds(ii,1):m_inds(ii,2);
    if icdiag == 1
        norm_parts(ii) = sum( mvec(inds).^2 .* Cmat(inds) );
    else
        norm_parts(ii) = mvec(inds)' * Cmat(inds,inds) * mvec(inds);  
    end
    disp(sprintf('%4i%16.4e%16.4e',ii,norm_parts(ii),sqrt(norm_parts(ii)) ));
end

if 1==1
    disp(sprintf('   1%16.4e',sum(norm_parts(1))));
    disp(sprintf(' 2-4%16.4e',sum(norm_parts(2:4))));
    disp(sprintf(' 1-4%16.4e',sum(norm_parts(1:4))));
end
    
%=========================================================

%   subroutine compute_norm_sq(imnorm, ievent_min, ievent_max, nevent, index_source, nmod, &
%         mvec, mvec_prior, cov_model, norm_sq_parts, norm_parts_weight)
% 
%     ! This computes the norm-squared of a model vector using the model covariance.
%     ! The dimensions of the input model vector are always the same, but this norm
%     ! takes into account the number of events used, which may be less than nevent.
%     ! UPDATE: If mvec_prior is present, then subtract from mvec: (m-mprior)^T Cm (m-mprior)
%     ! NOTE 1: mprior is only used is imnorm = 1
% 
%     integer, intent(in) :: imnorm, ievent_min, ievent_max, nevent, nmod
%     integer, dimension(NVAR_SOURCE, nevent), intent(in) ::index_source
%     double precision, dimension(nmod), intent(in) :: mvec, mvec_prior, cov_model
%     double precision, dimension(NVAR), intent(in), optional :: norm_parts_weight
% 
%     !double precision, intent(out) :: norm_total, norm_struct, norm_source
%     double precision, dimension(NVAR), intent(out) :: norm_sq_parts
% 
%     double precision, dimension(nmod) :: mtemp, ctemp
%     double precision, dimension(NVAR) :: npw
%     integer :: ievent, itemp1, itemp2, itemp3
% 
%     !----------
% 
%     norm_sq_parts(:) = 0.0 
%     ctemp(:) = 0.0
%     mtemp(:) = 0.0
%     npw(:) = 1.0
% 
%     ! NOTE 1: if taking the norm of a gradient, use the inverse covariance matrix
%     ! NOTE 2: if taking the norm of a model, use m - mprior
%     ! NOTE 3: norm_parts_weight is related to the norm-squared weights (not norm)
%     if (imnorm == 0) then
%         ctemp(:) = 1.0 / cov_model(:)
%         mtemp(:) = mvec(:)
%         if (present(norm_parts_weight)) npw(:) = 1.0 / norm_parts_weight(:)
%         !if (present(norm_parts_weight)) npw(:) = 1.0 / norm_parts_weight(:)**2
% 
%     elseif (imnorm == 1) then
%         ctemp(:) = cov_model(:)
%         mtemp(:) = mvec(:) - mvec_prior(:)
%         if (present(norm_parts_weight)) npw(:) = norm_parts_weight(:)
%         !if (present(norm_parts_weight)) npw(:) = norm_parts_weight(:)**2
% 
%     else
%         stop 'imnorm must = 0 or 1'
%     endif
% 
%     ! structure part of the norm -- BETA only
%     ! NOTE: division corresponds to inversion of a diagonal covariance matrix
%     norm_sq_parts(1) = sum( mtemp(1 : NLOCAL)**2 / ctemp(1 : NLOCAL) )
% 
%     ! source part of the norm -- only events that you are inverting for
%     do ievent = ievent_min, ievent_max
%        itemp1 = NLOCAL + index_source(1,ievent)
%        itemp2 = NLOCAL + index_source(2,ievent)
%        itemp3 = NLOCAL + index_source(3,ievent)
% 
%        norm_sq_parts(2) = norm_sq_parts(2) + mtemp(itemp1)**2 / ctemp(itemp1)
%        norm_sq_parts(3) = norm_sq_parts(3) + mtemp(itemp2)**2 / ctemp(itemp2)
%        norm_sq_parts(4) = norm_sq_parts(4) + mtemp(itemp3)**2 / ctemp(itemp3)
%     enddo
% 
%     !norm_struct = norm_sq_parts(1)
%     !norm_source = sum( norm_sq_parts(2:4) )
%     !norm_total  = sum( norm_sq_parts(1:4) )
% 
%     ! weight each part of the norm
%     norm_sq_parts(:) = norm_sq_parts(:) * npw(:)
% 
%   end subroutine compute_norm_sq
%   
%   