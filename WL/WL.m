clear
load ../ftn58sparse.mat
% ftn58sparse.dd = [ftn58sparse.dd zeros(size(ftn58sparse.dd,1),1)];
%% input parameter
%%-- Inital Setup --%%
ind.nkx = 300;  % the main direction of position operator
ind.nky = 300;
ind.nkz = 2;    %2
ind.nVB = 2; % # of occupy states (VB)
ind.Gvec = [1 0 0]*2*pi; % G-vecotr used for periodic gauge
ind.wsoc = 1;

nba = [ftn58sparse.Ainfo.Norb];
rr = reshape([ftn58sparse.Ainfo.Position],3,[])';
R=[];
for ir=1:length(nba)
    R=cat(1,R,repmat(rr(ir,:),[nba(ir) 1]));
end
ind.R=R;
nkz = ind.nkz;
nky = ind.nky;
%%-- set k-mesh --%%
ind.kx = linspace(-1*pi,1*pi,ind.nkx);
ind.kx = ind.kx(1:ind.nkx);
ind.ky = linspace(-1*pi,1*pi,ind.nky);
ind.kz = [0 pi]; %pi
%ind.kz = [0.2 0.4];
%%-- Load data --%%
ind.Norb = ftn58sparse.norb;
ind.ii   = ftn58sparse.ij(:,1);
ind.jj   = ftn58sparse.ij(:,2);
ind.tt   = ftn58sparse.tt;
ind.dd   = ftn58sparse.dd;
%%
phase = zeros(ind.nVB,nky,nkz);
CkWiLp = zeros(nky,nkz);
for ikz= 1:nkz
    for iky=1:nky
        fprintf('tot by=%d , count=%d \n',nky,iky);
        [fphase,fCkWiLp]=BPhase3D(ind,iky,ikz);
        phase(:,iky,ikz)=fphase;
        CkWiLp(iky,ikz)=fCkWiLp;
    end
end
%%
axisX = ind.ky/pi;
CkWiLp= CkWiLp/pi;
phase= phase/pi;
save BP.mat CkWiLp axisX phase ind
for ikz = 1:nkz
    C = 0;
    for n = 1:ind.nVB
        th = squeeze(phase(n,:,ikz)) * pi;   % 取第 n 個佔據帶的 Wilson loop 相位 (rad)
        th_unw = unwrap(th);                 % 解開相位跳躍
        C = C + (th_unw(end)-th_unw(1)) / (2*pi);
    end
    C = round(C);  % 四捨五入到最近整數，避免數值誤差
    fprintf('Chern number at kz index %d = %d\n', ikz, C);
end
eval(['run plot_WL.m'])
