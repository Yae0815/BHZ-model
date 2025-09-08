 function [ftn58_brick,ftn58_hop]=brickhop(ftn58sparse,surface,hop_d)

ddtmp=[ftn58sparse.dd];
ijtmp=[ftn58sparse.ij];
tttmp=[ftn58sparse.tt];

hopinfo=ddtmp(:,surface);

if hop_d==1
    brick_ind=find(hopinfo<1 & hopinfo>=0);
    hop_ind=find(ddtmp(:,surface)>0);
elseif hop_d ==-1
    brick_ind=find(hopinfo>-1 & hopinfo<=0);
    hop_ind=find(ddtmp(:,surface)<0);
end
ftn58_brick = ftn58sparse;
ftn58_brick.dd = ddtmp(brick_ind,:);
ftn58_brick.ij = ijtmp(brick_ind,:);
ftn58_brick.tt = tttmp(brick_ind,:);

ftn58_hop = ftn58sparse;
ftn58_hop.dd = ddtmp(hop_ind,:);
ftn58_hop.ij = ijtmp(hop_ind,:);
ftn58_hop.tt = tttmp(hop_ind,:);

end