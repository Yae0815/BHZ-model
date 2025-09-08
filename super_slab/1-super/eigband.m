function [eigenvalues, eigvec]=eigband(ftn58,kvector)

%global norb Hfix icross ftn58 subindex

norb=ftn58(1,1);
nbond=ftn58(1,2);

iii=ftn58(2:end,2);
%onsites=readcoordi('coordi.dat');
% Hfix=zeros(norb,norb);
% 
% for ib=1:nbond
%     if ftn58(ib+1,5)==0
%         if ftn58(ib+1,6)==0
%             if ftn58(ib+1,2)==ftn58(ib+1,3)
%             Hfix(ftn58(ib+1,2),ftn58(ib+1,3))=Hfix(ftn58(ib+1,2),ftn58(ib+1,3))+ftn58(ib+1,4)/2;
%         else
%             Hfix(ftn58(ib+1,2),ftn58(ib+1,3))=Hfix(ftn58(ib+1,2),ftn58(ib+1,3))+ftn58(ib+1,4);
%         end
%         end
%     end
% end
% Hfix=Hfix+Hfix';
% 
% icross=find(or(ftn58(2:end,5),ftn58(2:end,6)));
% 
% H=Hfix;
%  for ibb=icross'+1
%      ttemp=ftn58(ibb,4)*exp(i*(kvector(1)*ftn58(ibb,5)+kvector(2)*ftn58(ibb,6)));
%      H(ftn58(ibb,2),ftn58(ibb,3))=H(ftn58(ibb,2),ftn58(ibb,3))+ttemp;
%      H(ftn58(ibb,3),ftn58(ibb,2))=H(ftn58(ibb,3),ftn58(ibb,2))+conj(ttemp);
%  end
H=zeros(norb,norb);
 for ib=1:nbond
     ibb=ib+1;
     ttemp=ftn58(ibb,4)*exp(i*(kvector(1)*ftn58(ibb,5)+kvector(2)*ftn58(ibb,6)+kvector(3)*ftn58(ibb,7)));
     H(ftn58(ibb,2),ftn58(ibb,3))=H(ftn58(ibb,2),ftn58(ibb,3))+ttemp;
%     if ftn58(ibb,2)~=ftn58(ibb,3)
%     H(ftn58(ibb,3),ftn58(ibb,2))=H(ftn58(ibb,3),ftn58(ibb,2))+conj(ttemp);
%     end
 end

H=(H+H')/2;

%H=H(subindex,subindex);
%eigenvalues=eig(H);
[eigvec, eigenvalues]=eig(H);
eigenvalues=real(diag(eigenvalues));
return

