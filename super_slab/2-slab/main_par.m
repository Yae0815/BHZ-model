%function [ output_args ] = main1( input_args )
%MAIN1 Summary of this function goes here
%  Detailed explanation goes here
clear all

tic;

load Sftn58sparse.mat 
ftn58sparse=Sftn58sparse; 

%%---kpath---%%
nk = 100;
p1 = [linspace(0,0,nk)' linspace(0.5,0,nk)' linspace(0,0,nk)'];
p2 = [linspace(0,2/3,nk)' linspace(0,-1/3,nk)' linspace(0,0,nk)'];
kpts = [p1(1:nk,:);p2(1:nk,:)]*2*pi;

nk=length(kpts);
nb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
tt   = ftn58sparse.tt;
dd   = ftn58sparse.dd ;
orps = ftn58sparse.Orbitps(:,4:6);
Ek=zeros(nk,nb);
for ikk=1:nk
    fprintf('total k=%d now=%d\n',nk,ikk);
    kpoints=kpts(ikk,:);
    Hsparse = sparse(ii,jj,exp(1i*dd*kpoints').*tt,nb,nb);
    H0      = full(Hsparse);
    H0      = (H0 + H0')/2;
    [V,D]   = eig(H0);
    Etmp    = diag(D);
    Ek(ikk,:)=Etmp;
end

A=Ek;

toc;
Computation_time=num2str(toc); % computation time 

%save ftn59dope0.0cut0.0001.mat ftn58
%save Ek.mat A
%save time.mat Computation_time

%return

%save ftn58.mat ftn58
save Ek.mat A

return

plot(A,'r','LineWidth',2)
hold on
axis square


% Orbital
%dz2, dxz, dyz, dx2-y2, dxy
%pz, px, py

%on-site energy
%ii=find(ftn58(:,2)==ftn58(:,3) & ftn58(:,5)==0 & ftn58(:,6)==0 & ftn58(:,7)==0);ftn58(ii,2:7)

%hopping
%ii=find( ftn58(:,5)==0 & ftn58(:,6)==0 & ftn58(:,7)==0);ftn58(ii(2:end),2:7) 
