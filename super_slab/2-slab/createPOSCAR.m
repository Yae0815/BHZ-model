clear
load Sftn58sparse.mat
for i =1:Sftn58sparse.Nat
    tmp1=getfield(Sftn58sparse.Ainfo, {i},'Position');
    tmp2=getfield(Sftn58sparse.Ainfo, {i},'Atom');
    TT(i,1:3)= tmp1;
    if tmp2 == 'Mn'
        ttmp2 = 1;
    elseif tmp2 =='Bi'
        ttmp2 = 2;
    elseif tmp2 == 'Te'
        ttmp2 = 3;
    end
    A1(i)= ttmp2;
end
Ap = Sftn58sparse.BR3D;
[a b]=sort(A1);
numMn=length(find(a==1));
numTe=length(find(a==2));
numBi=length(find(a==3));
numatom = [numMn numTe numBi];
for j = 1:Sftn58sparse.Nat
    reT(j,1:3)= TT(b(j),:);
end

fid = fopen('creatPOSCAR.vasp','w+');
fprintf(fid,'creat-slab-MnBiTe \n');
fprintf(fid,'1.000 \n');
fprintf(fid,' %8.7f %8.7f %8.7f \n',Ap');
fprintf(fid,'Mn Te Bi \n');
%%
fprintf(fid,' %3.0f %3.0f %3.0f \n',numatom');
fprintf(fid,' Rec \n');
fprintf(fid,' %8.7f  %8.7f  %8.7f \n',reT');
fclose all;
