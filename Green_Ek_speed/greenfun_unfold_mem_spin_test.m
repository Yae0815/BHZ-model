%%% Green's function tech. %%%
%%%    ====   Ek   ====    %%%
%%% v1. Tay-Rong Chang     %%%
%%% v2. Yi-Hsin Tu         %%%
%%% v3. Yueh-Ting Yao      %%%

%%% (1) fold and unfold       %%%
%%% (2) semi-infinite SS      %%%
%%% (3) spin decompose        %%%
%%% (4) orbital decompose     %%%
%%% (5) choose surface atoms  %%%
clear 
% c = parcluster('local');
% c.NumWorkers = 16;
% parpool(c, c.NumWorkers);


load ../ftn58sparse.mat
% ftn58sparse = SPftn58sparse;
lat = ftn58sparse.BR.*ftn58sparse.abc;
rec = inv(lat)';

%---kpath---%
nk = 500;
p1 = [linspace(0.0,0.0,nk)' linspace(0.5,0.0,nk)' linspace(0,0,nk)'];
p2 = [linspace(0.0,0.0,nk)' linspace(0.0,-0.5,nk)' linspace(0,0,nk)'];
kpts = [p1(1:nk,:);p2(1:nk,:)]*2*pi;


%%
tic; % start clock
ticBytes(gcp);
%% ===input start===%
Ni=10;                          % # of interation (layer); layer=2^Ni
EF=0;                           % Fermi energy
iter=1;                         % iteration 1 or 0
wsoc=1;                         % 1 for w/ soc; 0 for w/o soc
w=EF+linspace(-1,1,1000);       % energy window %1000
nn=0.001;                       % energy broadening

surface=3; % x=1,y=2,z=3
hop_d=-1; % bottom surface = 1, upper surface = -1

%%% unfolding %%%
% Atom={'W1','W1','W2','W2','Se1','Se1','Se2','Se2','Se3','Se3','Se4','Se4'};
Atom={};
for ti = 1:ftn58sparse.Nat
    Atom = [Atom,num2str(ti)];
end

%%% orbital selection %%%
nb=ftn58sparse.norb;
oo=[1:nb];                  % orbital selection 

%%% plot surface %%%
%sur_orb = 2;                % # of surface orbital
sur_orb = 4   %測試調整
orbit = ftn58sparse.Orbitps;
if hop_d == -1
    tmp = sortrows(orbit,surface+3,'descend'); % upper surface = -1
elseif hop_d == 1
    tmp = sortrows(orbit,surface+3,'ascend'); % bottom surface = 1
end
orbit_plot = tmp(1:sur_orb,1);

% orbital_onsite = [ 31:36 43:45 ];
% orbital_onsite = [orbital_onsite orbital_onsite+nb/2];
% orbit = ftn58sparse.Orbitps;
% tmp = sortrows(orbit,surface+3,'descend'); % upper surface = -1
% tmp = sortrows(orbit,surface+3,'ascend'); % bottom surface = 1
% orbit_plot = tmp(1:54,1);
% orbital_onsite = tmp(1:24,1);
%---kpath---%
% nk = 100;
% p1 = [linspace(0,0,nk+1)' linspace(0.5,0,nk+1)' linspace(0,0,nk+1)']; % M G
% p2 = [linspace(0,2/3,nk)' linspace(0,-1/3,nk)' linspace(0,0,nk)'];    % G K
% kpts = [p1(1:nk,:);p2(1:nk,:)]*2*pi;
%---kmesh(CEC)---%
% LatM=[ftn58sparse.BR].*ftn58sparse.abc;
% recM=[cross(LatM(2,:),LatM(3,:));cross(LatM(3,:),LatM(1,:));cross(LatM(1,:),LatM(2,:))]/dot(LatM(1,:),cross(LatM(2,:),LatM(3,:)));
% recM=[recM(1,:)/norm(recM(1,:));recM(2,:)/norm(recM(2,:));recM(3,:)/norm(recM(3,:))];
% 
% 
% nk = [10 1 10];
% kkx = linspace(-0.01,0.01,nk(1));
% kky = linspace(0,0,nk(2));
% kkz = linspace(-0.2,0.2,nk(3));
% 
% [ky,kx,kz]=meshgrid(kky,kkx,kkz);
% kpts = [kx(:),ky(:),kz(:)]*inv(recM)*2*pi;
%% ===input end===%%
nk=length(kpts);
nb=ftn58sparse.norb;
Ek=zeros(nk,nb);
%---unfold-Matric---%
nba = [ftn58sparse.Ainfo.Norb];
rr = reshape([ftn58sparse.Ainfo.Position],3,[])';
M=unfoldope(Atom,nba,wsoc);
R=[];
for ir=1:length(nba)
    R=cat(1,R,repmat(rr(ir,:),[nba(ir) 1]));
end
%---------------------%
II=eye(nb);
AA=zeros(length(w),nk);
%---------------------% spin
No  = nb/2;          % 軌域數
I_o = eye(No);       % τ₀
sx  = kron([0 1; 1 0],  I_o);   % τ₀⊗σₓ
sy  = kron([0 -1i; 1i 0], I_o); % τ₀⊗σᵧ
sz  = kron([1  0; 0 -1],  I_o); % τ₀⊗σ_z

%---------------------%
Ab=AA;
Abx=AA;
Aby=AA;
Abz=AA;
Aux=AA;
Auy=AA;
Auz=AA;
As=AA;
Asfold=AA;
Asur=AA;
Afold=AA;
Aunfold=AA;
Asx=AA;
Asy=AA;
Asz=AA;
% As_part = zeros(length(w),nk);
%---------------------%
if iter == 0
    [HH_H_mat_4,HH_dd_hop]=reshape_ftn58(ftn58sparse);
    for ik=1:nk
        fprintf('total k = %d now = %d\n',nk,ik);
%         [Hk] = genes_mem(ftn58sparse,kpts(ik,:),wsoc,R,0);
        [Hk] = genes_mem_no_sparse(HH_H_mat_4,HH_dd_hop,kpts(ik,:),wsoc,R,0,nb);
        Afoldt=zeros(1,length(w));
        Aunt=zeros(1,length(w));
        Aunxt=zeros(1,length(w));
        Aunyt=zeros(1,length(w));
        Aunzt=zeros(1,length(w));
        Hw=Hk(:,:);
        for iw=1:length(w)
            runw=w(iw);
            E=(runw+1i*nn)*II;    % w in geen function
            GG=inv(E-Hw);

            Afoldt(1,iw)=-imag(trace(GG));   %% folding

            MGG = M*GG;
             MGbx=sx*MGG; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             MGby=sy*MGG; % comment if you have without SOC case %%
             MGbz=sz*MGG; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Aunt(1,iw)=-imag(trace(MGG));   % spectral weight of bulk
            Aunxt(1,iw)=-imag(trace(MGbx));
            Aunyt(1,iw)=-imag(trace(MGby));
            Aunzt(1,iw)=-imag(trace(MGbz));

            tmp_part=-imag(diag(M*GG));
        end
        Afold(:,ik)=Afoldt;
        Aunfold(:,ik)=Aunt;
        Aux(:,ik)=Aunxt;
        Auy(:,ik)=Aunyt;
        Auz(:,ik)=Aunzt;
    end
elseif iter == 1
    [HH,TT]=brickhop(ftn58sparse,surface,hop_d);
    [HH_H_mat_4,HH_dd_hop]=reshape_ftn58(HH);
    [TT_H_mat_4,TT_dd_hop]=reshape_ftn58(TT);
%     [Hsur,Tsur]=brickhop(ftn58sparse,surface,hop_d);
% clear ftn58sparse
    parfor ik=1:nk
        fprintf('total k = %d now = %d\n',nk,ik);
%         [Hk] = genes_mem(HH,kpts(ik,:),wsoc,R,0);
%         [Tk] = genes_mem(TT,kpts(ik,:),wsoc,R,1);
        [Hk] = genes_mem_no_sparse(HH_H_mat_4,HH_dd_hop,kpts(ik,:),wsoc,R,0,nb)
        [Tk] = genes_mem_no_sparse(TT_H_mat_4,TT_dd_hop,kpts(ik,:),wsoc,R,1,nb)
%         [Hs0] = genes_mem(Hsur,kpts(ik,:),wsoc,R,0);
        Hs0 = Hk;
%         [Tk0] = genes_mem(Tsur,kpts(ik,:),wsoc,R,1);
        Tk0 = Tk;
        Abt=zeros(1,length(w));
        Abxt=zeros(1,length(w));
        Abyt=zeros(1,length(w));
        Abzt=zeros(1,length(w));
        Ast=zeros(1,length(w));
        Asfoldt=zeros(1,length(w));
        Asxt=zeros(1,length(w));
        Asyt=zeros(1,length(w));
        Aszt=zeros(1,length(w));
        Asurt=zeros(1,length(w));
        Ast_part=zeros(1,length(w));
        H=Hk(:,:);
        T=Tk(:,:);
        H0=Hs0(:,:);
        for iw=1:length(w)
            runw=w(iw);
            E=(runw+1i*nn)*II;    % w in geen function  %same
%             [Hs,Hb]=renormH(E,H,T,Ni);  % PRB 31, 5166 eq.10
            [Hs,Hb]=renormH_speed(E,H,T,Ni);  % PRB 31, 5166 eq.10
            Gb=inv(E-Hb);  % bulk  % PRB 31, 5166 eq.5
            Gs=inv(E-Hs);  % surface
            Gs0=inv(E-H0-T*Gs*T');    % cover surface % PRB 31, 5166 eq.13

            MGs0=M*Gs0;
             MGsx=sx*MGs0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             MGsy=sy*MGs0; % comment if you have without SOC case %%
             MGsz=sz*MGs0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            MGb=M*Gb;
             MGbx=sx*MGb; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             MGby=sy*MGb; % comment if you have without SOC case %%
             MGbz=sz*MGb; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Abt(1,iw)=-imag(trace(MGb));   % spectral weight of bulk
            Abxt(1,iw)=-imag(trace(MGbx));
            Abyt(1,iw)=-imag(trace(MGby));
            Abzt(1,iw)=-imag(trace(MGbz));

            Ast(1,iw)=-imag(trace(MGs0));
			Asfoldt(1,iw)=-imag(trace(Gs0));
             Asxt(1,iw)=-imag(trace(MGsx));
             Asyt(1,iw)=-imag(trace(MGsy));
             Aszt(1,iw)=-imag(trace(MGsz));
            tmp_part=-imag(diag(MGs0));
%             tmp_part=-imag(diag(MGs0));
            Asurt(1,iw)=-imag(trace(MGs0(oo,oo)));
            Ast_part(1,iw) = sum(tmp_part(orbit_plot));
        end
        Ab(:,ik)=Abt;
        Abx(:,ik)=Abxt;
        Aby(:,ik)=Abyt;
        Abz(:,ik)=Abzt;
        
        As(:,ik)=Ast;
        Absfold(:,ik)=Asfoldt;
        Asx(:,ik)=Asxt;
        Asy(:,ik)=Asyt;
        Asz(:,ik)=Aszt;
        Asur(:,ik)=Asurt;
        As_part(:,ik)=Ast_part;
    end
end
% ----------------------------------- %
save ek.mat w nk EF A* nb -v7.3
% save ek_-1_3_no_unfold_-X_G_X.mat w nk EF A* orbit tmp orbit_plot orbital_onsite nb
toc; % stop clock
tocBytes(gcp)
