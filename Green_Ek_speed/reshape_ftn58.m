 function [H_mat_4,dd_hop] = reshape_ftn58(ftn58sparse)

nb = ftn58sparse.norb;
ij = [ftn58sparse.ij];
ii   = ij(:,1);
jj   = ij(:,2);
tt   = ftn58sparse.tt;
dd   = ftn58sparse.dd;

dd_hop = unique(dd,'rows'); 
dd_hop_2 = [dd_hop dd_hop(:,1)*10000+dd_hop(:,2)*100+dd_hop(:,3)*1];

H_mat_2 = cell(1,size(dd_hop,1));
H_zero = zeros(nb,nb);
for i = 1:size(dd_hop,1)
    H_mat_2{1,i} = H_zero;  
end
dd_index = dd(:,1)*10000+dd(:,2)*100+dd(:,3)*1;
[~,index] = ismember(dd_index,dd_hop_2(:,4)); 
for i = 1:size(dd,1)
%     dd_index = dd(i,1)*10000+dd(i,2)*100+dd(i,3)*1; % takes long
%     index = find(dd_index==dd_hop_2(:,4)); % takes long
%     [~,index] = ismember(dd_hop_2(:,1:3),dd(i,:),'rows'); % takes very long
%     H_mat_2{1,index}(ii(i),jj(i)) = H_mat_2{1,index}(ii(i),jj(i)) + tt(i);
    H_mat_2{1,index(i)}(ii(i),jj(i)) = H_mat_2{1,index(i)}(ii(i),jj(i)) + tt(i);
end
H_mat_3 = zeros(nb,nb,size(dd_hop,1));
for i = 1:size(dd_hop,1)
    H_mat_3(:,:,i) = H_mat_2{1,i};
end
H_mat_4 = reshape(H_mat_3,nb^2,size(dd_hop,1));

end