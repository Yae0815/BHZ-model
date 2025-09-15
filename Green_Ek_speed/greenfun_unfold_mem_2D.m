%% ==== Constant-Energy Cut (E = EF) on 2D k-mesh ====
clear; clc;

% ---------- Load model ----------
load ../ftn58sparse.mat

% ---------- User controls ----------
EF      = 0.0;        % Fermi level (eV)
nn      = 1e-4;       % broadening (eV)
wsoc    = 1;          % 1 with SOC, 0 without
surface = 3;          % stacking (cleavage) axis: x=1, y=2, z=3
hop_d   = 1;          % surface side: bottom=+1, top=-1
Ni      = 10;          % decimation depth (2^Ni layers)

% 2D k-mesh (fractional along reciprocal primitive axes)
Nkx = 201; 
Nky = 201;
kr1 = [-0.5, 0.5];    % range along the first FREE axis (in 2π units)
kr2 = [-0.5, 0.5];    % range along the second FREE axis (in 2π units)
kfix_frac = 0.03;     % fixed fraction along the STACKING axis

% what to plot: {'As','Asx','Asy','Asz','As_part','Ab'}
target = 'As';
sur_orb = 4;          % for As_part: how many topmost orbitals to sum

% ---------- Dimensions / projectors ----------
nb   = ftn58sparse.norb;
Atom = arrayfun(@(i) num2str(i), 1:ftn58sparse.Nat, 'uni', 0);
nba  = [ftn58sparse.Ainfo.Norb];
rr   = reshape([ftn58sparse.Ainfo.Position], 3, [])';  % Nat x 3 (fractional direct)

% Unfolding projector
M = unfoldope(Atom, nba, wsoc);

% Build R: repeat each atom's fractional position for its orbitals
R = [];
for ir = 1:length(nba)
    R = cat(1, R, repmat(rr(ir,:), [nba(ir) 1]));  % nb x 3
end

II = eye(nb);

% Spin matrices (only meaningful when wsoc==1 and nb is even)
if mod(nb,2)==0 && wsoc==1
    dnup = eye(nb/2); updn = dnup;
    sx = [zeros(nb/2), dnup; updn, zeros(nb/2)];
    dnup = -1i*eye(nb/2); updn = -dnup;
    sy = [zeros(nb/2), dnup; updn, zeros(nb/2)];
    upup = eye(nb/2); dndn = -upup;
    sz = [upup, zeros(nb/2); zeros(nb/2), dndn];
else
    sx = zeros(nb); sy = zeros(nb); sz = zeros(nb);
end

% ---------- Select top-layer orbitals for As_part ----------
orbit = ftn58sparse.Orbitps;   % columns include fractional x y z at (col 4~6)
if hop_d == -1
    tmp = sortrows(orbit, surface+3, 'descend'); % top surface
else
    tmp = sortrows(orbit, surface+3, 'ascend');  % bottom surface
end
sur_orb = min(sur_orb, size(tmp,1));
orbit_plot = tmp(1:sur_orb,1);

% ---------- Build 2D k-mesh (auto based on surface axis) ----------
free = setdiff(1:3, surface);  % the two in-plane reciprocal axes
[k1g, k2g] = meshgrid(linspace(kr1(1), kr1(2), Nkx), ...
                      linspace(kr2(1), kr2(2), Nky));
Ntot = numel(k1g);

% fractional k in units of reciprocal basis (b1,b2,b3); multiply by 2π later
kfrac = zeros(Ntot,3);
kfrac(:, free(1)) = k1g(:);
kfrac(:, free(2)) = k2g(:);
kfrac(:, surface) = kfix_frac;

% match your original convention: pass kpts in "fractional * 2*pi"
kpts = kfrac * (2*pi);   % Ntot x 3

% ---------- Semi-infinite surface GF blocks ----------
[HH,TT] = brickhop(ftn58sparse, surface, hop_d);
[HH_H_mat_4,HH_dd_hop] = reshape_ftn58(HH);
[TT_H_mat_4,TT_dd_hop] = reshape_ftn58(TT);

% ---------- Allocate ----------
As_E  = zeros(Ntot,1);    % surface total spectral weight at EF
Asx_E = zeros(Ntot,1);    % spin-resolved
Asy_E = zeros(Ntot,1);
Asz_E = zeros(Ntot,1);
AsP_E = zeros(Ntot,1);    % partial (top-layer orbitals)
Ab_E  = zeros(Ntot,1);    % bulk (optional)

E = (EF + 1i*nn) * II;

% ---------- (Optional) quick sanity check on one k before parfor ----------
% ik_test = 1;
% HkT = genes_mem_no_sparse(HH_H_mat_4,HH_dd_hop,kpts(ik_test,:),wsoc,R,0,nb);
% TkT = genes_mem_no_sparse(TT_H_mat_4,TT_dd_hop,kpts(ik_test,:),wsoc,R,1,nb);

% ---------- Parallel over k ----------
% 若需要平行，解除以下註解
% c = gcp('nocreate'); if isempty(c), parpool('local'); end

parfor ik = 1:Ntot
    % k-dependent Hamiltonians (IMPORTANT: pass R)
    Hk = genes_mem_no_sparse(HH_H_mat_4,HH_dd_hop,kpts(ik,:),wsoc,R,0,nb);
    Tk = genes_mem_no_sparse(TT_H_mat_4,TT_dd_hop,kpts(ik,:),wsoc,R,1,nb);

    % decimation (PRB 31, 5166)
    [Hs,Hb] = renormH_speed(E,Hk,Tk,Ni);

    % surface Green's function with capping (Eq.13)
    Gs  = inv(E - Hs);
    Gs0 = inv(E - Hk - Tk*Gs*Tk');   % cover surface

    % bulk Green's function (optional)
    Gb = inv(E - Hb);

    % unfolding projections
    MGs0 = M*Gs0;     
    MGb  = M*Gb;

    % surface totals & spin components
    As_E(ik)  = -imag(trace(MGs0));
    Asx_E(ik) = -imag(trace(sx*MGs0));
    Asy_E(ik) = -imag(trace(sy*MGs0));
    Asz_E(ik) = -imag(trace(sz*MGs0));

    % bulk (optional)
    Ab_E(ik) = -imag(trace(MGb));

    % partial: sum of local diagonal on chosen orbitals
    diag_surf = -imag(diag(MGs0));
    AsP_E(ik) = sum(diag_surf(orbit_plot));
end

% ---------- Choose target map ----------
switch lower(target)
    case 'as',      val = As_E;
    case 'asx',     val = Asx_E;
    case 'asy',     val = Asy_E;
    case 'asz',     val = Asz_E;
    case 'as_part', val = AsP_E;
    case 'ab',      val = Ab_E;
    otherwise,      val = As_E;
end

Map = reshape(val, size(k1g)) / pi;   % (1/pi) normalization as usual

%% ---------- Plot (pcolor + shading interp, viridis, LaTeX labels) ----------
% 參數：你可以改這兩個來控制外觀
cax = [0 10];           % 顏色範圍 caxis，依資料調整
tick_step = 0.25;       % k 軸刻度間距（以 2π 分數單位）

% 1) 用 pcolor + shading interp 畫（跟 ek.mat 風格一致）
figure('Color','w');
hold on;
% pcolor 需要格點；我們直接用中心格點即可（會少掉最後一行／列是正常現象）
hp = pcolor(k1g, k2g, Map);
set(hp,'EdgeColor','none');    % 等效 shading interp 的平滑感
shading interp;

% colormap / 色條
try
    colormap(viridis);  % 你資料夾有 viridis.m 就會用到它
catch
    colormap(parula);
end
cb = colorbar;
% 針對不同 target 給漂亮一點的標籤
switch lower(target)
    case 'as',      cb.Label.String = '$A_s/\pi$ [1/eV]';
    case 'asx',     cb.Label.String = '$A_{s,x}/\pi$ [1/eV]';
    case 'asy',     cb.Label.String = '$A_{s,y}/\pi$ [1/eV]';
    case 'asz',     cb.Label.String = '$A_{s,z}/\pi$ [1/eV]';
    case 'as_part', cb.Label.String = '$A_{s}^{(\mathrm{top})}/\pi$ [1/eV]';
    case 'ab',      cb.Label.String = '$A_b/\pi$ [1/eV]';
    otherwise,      cb.Label.String = sprintf('%s / \\pi  [1/eV]', upper(target));
end
cb.Label.Interpreter = 'latex';
cb.TickDirection = 'out';
cb.FontSize = 14;

% 座標與比例
axis image; axis xy;
xlim(kr1); ylim(kr2);
caxis(cax);

% 2) 讓軸標顯示你想要的「x / y / z」
%    free = setdiff(1:3, surface) 已在主程式前面算好了
axes_lbl_xyz = {'x','y','z'};    % 想用 b_1/b_2/b_3 的話換下一行
% axes_lbl_b   = {'b_1','b_2','b_3'};

% 用 LaTeX interpreter 才能正常顯示下標 k_x, k_y
xlabel(sprintf('$k_{%s}$ (fraction of $2\\pi$)', axes_lbl_xyz{free(1)}), 'Interpreter','latex');
ylabel(sprintf('$k_{%s}$ (fraction of $2\\pi$)', axes_lbl_xyz{free(2)}), 'Interpreter','latex');

% 標題（只放 EF 與固定軸）
title(sprintf('$E=E_F=%.3f\\,\\mathrm{eV}$;  fixed $k_{%s}=%.3f$', ...
      EF, axes_lbl_xyz{surface}, kfix_frac), 'Interpreter','latex');

% 刻度 & 風格（比照 ek.mat）
ax = gca;
ax.TickDir    = 'in';
ax.FontSize   = 18;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.TickLabelInterpreter = 'latex';

% 設定等距刻度
xticks(kr1(1):tick_step:kr1(2));
yticks(kr2(1):tick_step:kr2(2));
box on;
hold off;

% ---------- Save ----------
save cec_EF_2D.mat EF surface hop_d kfix_frac free kr1 kr2 Nkx Nky ...
     target Map k1g k2g As_E Asx_E Asy_E Asz_E AsP_E Ab_E
