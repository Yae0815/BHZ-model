function [A,vecnorm]=bandplot(ftn58, kpoints)

%global norb Hfix icross ftn58 subindex


norb=ftn58(1);
kmn=size(kpoints,1);
A=zeros(kmn*norb,2);
eigvec=zeros(norb,norb,kmn);
for ik=1:kmn
    A(ik*norb-norb+1:ik*norb,1)=ik;
    fprintf('total k=%d now=%d\n',kmn,ik);
%     kx=klist(ik,1)/klist(ik,4)*2*pi;
%     ky=klist(ik,2)/klist(ik,4)*2*pi;
%     kz=klist(ik,3)/klist(ik,4)*2*pi;
    A(ik*norb-norb+1:ik*norb,2)=eigband(ftn58,kpoints(ik,:)*2*pi);
    [eigenvalues, eigvectemp]=eigband(ftn58,kpoints(ik,:)*2*pi);
    eigvec(:,:,ik)=eigvectemp;
end
    vecnorm=abs(eigvec).^2;
return
