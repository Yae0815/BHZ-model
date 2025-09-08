clear all

load super_ftn58sparse.mat;
% load ftn58sparse.mat;
ib1 = SPftn58sparse.ij(:,1);
ib2 = SPftn58sparse.ij(:,2);
% ib1 = ftn58sparse.ij(:,1);
% ib2 = ftn58sparse.ij(:,2);
%%
temp1 = find(ib1==1);
aa=ib2(temp1);
for i=1:180
temp(i,1)=size(find(aa==i),1);    
% temp2 = find(aa==85);
% temp3 = find(aa==90);
end
% bb=aa(temp2);
% cc=ib1(temp2);
% dd=ib1(temp1);

