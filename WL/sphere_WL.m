clear
load ../ftn58sparse.mat

%% ---- input / mesh ----
ind.nVB   = 2;
ind.Norb  = ftn58sparse.norb;
ind.ii    = ftn58sparse.ij(:,1);
ind.jj    = ftn58sparse.ij(:,2);
ind.tt    = ftn58sparse.tt;
ind.dd    = ftn58sparse.dd;

% sphere center & radius (把小球中心放在疑似 Weyl 點附近)
k0 = [-0.0456*2*pi, 0.1222*2*pi, 0];           % 你可以換成找到的 Weyl 位置
%k0 = [0, 0, 0];
Rk = 0.05*pi;             % 半徑要小到只包一顆 Weyl，且遠離其他簡并
Ntheta = 300; Nphi = 300;
theta = linspace(1e-3, pi-1e-3, Ntheta);   % 避免極點
phi   = linspace(0, 2*pi, Nphi+1); phi(end)=[];  % 週期

% 預存佔據子空間
Vocc = cell(Ntheta, Nphi);

% 掃球面所有點：求佔據本徵向量
for it = 1:Ntheta
    st = sin(theta(it)); ct = cos(theta(it));
    for ip = 1:Nphi
        cp = cos(phi(ip)); sp = sin(phi(ip));
        k = k0 + Rk*[st*cp, st*sp, ct];

        % H(k)
        Hs = sparse(ind.ii, ind.jj, exp(1i*(ind.dd*k.')).*ind.tt, ind.Norb, ind.Norb);
        H  = (Hs + Hs')/2;                 % 保 Hermitian

        [V, D] = eig(full(H), 'vector');
        [~, idx] = sort(real(D), 'ascend'); % 能量升冪
        W = V(:, idx(1:ind.nVB));
        Vocc{it, ip} = orth(W);             % 正交化，數值穩定
    end
end

%% ---- (A) Chern number on S^2 via FHS ----
% link variables (det-overlap) on theta,phi directions
Utheta = zeros(Ntheta-1, Nphi);
Uphi   = zeros(Ntheta,   Nphi);

% ---- theta-links ----
for it = 1:Ntheta-1
  for ip = 1:Nphi
    M = Vocc{it,ip}' * Vocc{it+1,ip};
    [U,S,V] = svd(M,'econ');
    L = U*V';                   % polar unitary of M
    Utheta(it,ip) = det(L);     % 已在單位圓上
    % 可選：退避小奇異值
    if min(diag(S)) < 1e-6, Utheta(it,ip) = Utheta(it,ip)/abs(Utheta(it,ip)); end
  end
end

% ---- phi-links ----
for it = 1:Ntheta
  for ip = 1:Nphi
    jp = ip+1; if jp > Nphi, jp = 1; end
    M = Vocc{it,ip}' * Vocc{it,jp};
    [U,S,V] = svd(M,'econ');
    L = U*V';
    Uphi(it,ip) = det(L);
    if min(diag(S)) < 1e-6, Uphi(it,ip) = Uphi(it,ip)/abs(Uphi(it,ip)); end
  end
end

% plaquette curvature phases and Chern
Fsum = 0.0;
for it = 1:Ntheta-1
    for ip = 1:Nphi
        jp = ip + 1; if jp > Nphi, jp = 1; end
        plaq = Utheta(it,ip) * Uphi(it+1,ip) * ...
               conj(Utheta(it,jp)) * conj(Uphi(it,ip));
        Fsum = Fsum + angle(plaq);
    end
end
Chern_S2 = round(Fsum/(2*pi));   % 最接近整數
fprintf('[S2-FHS] Chern (chiral charge) = %d\n', Chern_S2);

%% ---- (B) Wilson loop along phi at fixed theta ----
% 看相位流（等價於計 Weyl 單極子穿通數的纏繞）
nW = Nphi; phases = zeros(Ntheta, ind.nVB);
for it = 1:Ntheta
    % 建 φ-圈 Wilson 乘積
    F = eye(ind.nVB);
    for ip = 1:Nphi
        jp = ip + 1; if jp > Nphi, jp = 1; end
        M = Vocc{it,ip}' * Vocc{it,jp};
        % 極分解把非酉性去掉（比直接 det 更穩一點）
        [U,~,Vh] = svd(M);  L = U*Vh';
        F = F * L;
    end
    w = eig(F);
    phases(it, :) = sort(mod(angle(w), 2*pi)); % 0~2π
end

% 你可以把 phases(it,:) 對 theta 畫線，看是否有穿越（纏繞次數=手性）

save BP_S2.mat Chern_S2 phases theta phi k0 Rk

eval('run plot_sphere_WL.m');