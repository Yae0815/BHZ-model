function [phase,CkWiLp]=BPhase3D(ind,iky,ikz)
CkWiLp_tmp=1;
px_tmp=eye(ind.nVB,ind.nVB);
nkx = ind.nkx;
for ikx=1:nkx
    kpoints = [ind.kx(ikx) ind.ky(iky) ind.kz(ikz)]; %integral:x, along:y, specific:z
    %kpoints = [ind.kz(ikz) ind.kx(ikx) ind.ky(iky)]; %integral:z, along:x, specific:y
    %kpoints = [ind.ky(iky) ind.kz(ikz) ind.kx(ikx)]; %integral:y, along:z, specific:x
    U=kron(eye(ind.wsoc+1),diag(exp(-1i*ind.R*kpoints')));
    %%% Prepare F(N-1,0) %%%
    if ikx == 1
        Hsparse = sparse(ind.ii,ind.jj,exp(1i*ind.dd*kpoints').*ind.tt,ind.Norb,ind.Norb);
        HH = full(Hsparse);
        HH      = U*HH*U';
        HH = (HH+HH')/2;
        [V,E] = eig(HH);
        evec0 = V(:,1:ind.nVB);
        eveco = V(:,1:ind.nVB);
    else %%% ikx = 2:kx
        Hsparse = sparse(ind.ii,ind.jj,exp(1i*ind.dd*kpoints').*ind.tt,ind.Norb,ind.Norb);
        HH = full(Hsparse);
        HH      = U*HH*U';
        HH = (HH+HH')/2;
        [V,E] = eig(HH);
        evecn = V(:,1:ind.nVB);
        Ft = eveco'*evecn;
        px_tmp = px_tmp*Ft;
        CkWiLp_tmp= CkWiLp_tmp*det(Ft);
        eveco = evecn;
    end
end
%%% Prepare F(N,0)%%%
evecG = kron(eye(ind.wsoc+1),diag(exp(-1i*ind.R*ind.Gvec')))*evec0;  %% u_N periodic gauge (next BZ)
F0 = eveco'*evecG; %% F0 = u_0'*u_N
px_tmp=px_tmp*F0;
phase = angle(eig(px_tmp));
CkWiLp_tmp = CkWiLp_tmp*det(F0);
CkWiLp = angle(CkWiLp_tmp);