clear all

load Sftn58sparse.mat
ftn58sparse=Sftn58sparse;
ftn58sparse.dd=[Sftn58sparse.dd zeros(length(Sftn58sparse.dd(:,1)),1)];

ftn58=zeros(length(ftn58sparse.ij)+1,7);
ftn58(1,1)=ftn58sparse.norb;
ftn58(1,2)=length(ftn58sparse.ij);
ftn58(2:end,1)=[1:length(ftn58sparse.ij)]';
ftn58(2:end,2:3)=real(ftn58sparse.ij);
ftn58(2:end,4)=ftn58sparse.tt;
ftn58(2:end,5:7)=real(ftn58sparse.dd);

save ftn58.mat ftn58 -v7.3
