function [Hk]=genes_mem(ftn58sparse,kpts,wsoc,R,isT)
% nk=length(kpts);
nb=ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
tt   = ftn58sparse.tt;
dd   = ftn58sparse.dd;
% parfor ikk=1:nk
% for ikk=1:nk
% fprintf('total k=%d now=%d\n',nk,ikk);
kcc=kpts;
U=kron(eye(wsoc+1),diag(exp(-1i*R*kcc')));
Hsparse = sparse(ii,jj,exp(1i*dd*kcc').*tt,nb,nb);
H0      = full(Hsparse);
H0      = U*H0*U';
if isT==0
H0      = (H0 + H0')/2;
end
Hk(:,:) = H0;
% end
end
