function scan_DO_ftn58_consistent()
% BHZ lattice (spin ⊗ orbital). Only DO*tau_x (TRS kept, P broken).
% >>> ftn58-consistent mass sign <<<
% M_lat(k) = M + 2B1(cos kz - 1) + 2B2x(cos kx - 1) + 2B2y(cos ky - 1)
% Nodes satisfy on kz-plane:  M_lat=0  and  DO^2 = (A2x sin kx)^2 + (A2y sin ky)^2.
% We solve along M=0 contour (continuous, no eig, no k·p).
%
% Author: you+ChatGPT (2025-09-02)

%% ===== USER PARAMS (match ftn58!) =====
A1  = 0.2;
A2x = 0.2;    A2y = 0.2;               % 可維持等向；離散 Weyl 建議稍微各向異性 via B2x≠B2y
B1  = 0.3;
B2x = 0.3000;                           % 對應 ftn58 裡 cos kx 的係數
B2y = 0.3150*1.05;                      % 對應 ftn58 裡 cos ky 的係數（稍大 5%）
M    = -0.030;                          % 注意：這裡用 ftn58 的 M（非「M0 +⋯」）
C = 0; D1 = 0; D2 = 0;                  %#ok<NASGU>

KZ_PLANE = 0;                           % 设 0 或 pi 來選 kz 平面
PAD      = 0.01;                        % DO 視窗外擴 %
N_SCAN   = 4000;                        % 沿等值線的一維掃描密度
DO_MANUAL= [];                          % 若指定，直接用這個 DO 不算視窗

% 小盒子 3D 驗證
KZ_RANGE = 0.15*pi; NK_LOC = 25;

%% ===== 定義 kz 平面上的連續函數（ftn58 號誌） =====
cz = cos(KZ_PLANE);                     % kz=0 → cz=1；kz=π → cz=-1
% M=0 → cos ky = cy(kx)
Cy_from_kx = @(kx) 1 - ((-M) - 2*B1*(1-cz) - 2*B2x*(1-cos(kx))) / (2*B2y);
valid_ky   = @(kx) abs(Cy_from_kx(kx))<=1+1e-14;
ky_from_kx = @(kx) acos( clamp(Cy_from_kx(kx), -1, 1) );
R_of_kx    = @(kx) hypot( A2x*sin(kx), A2y*sin(ky_from_kx(kx)) ); % r(k) on M=0

%% ===== 找出 M=0 的有效 kx 區段 =====
ks = linspace(0, pi, N_SCAN+1);
mask = arrayfun(valid_ky, ks);
segs = contiguous_segments(ks(mask));
if isempty(segs)
    error('這組 (M,B1,B2x,B2y) 在 kz=%.2f 平面沒有 M=0 等值線；換平面或調 M。', KZ_PLANE);
end

%% ===== 連續的 DO 視窗（沿等值線的 min/max r） =====
r_min = +inf; r_max = -inf;
for s = 1:size(segs,1)
    a = segs(s,1); b = segs(s,2);
    rmin_s = fminbnd(@(x) R_of_kx(x), a, b);
    rmax_s = fminbnd(@(x) -R_of_kx(x), a, b);
    r_min = min(r_min, R_of_kx(rmin_s));
    r_max = max(r_max, R_of_kx(rmax_s));
end
if isempty(DO_MANUAL)
    DO_lo = r_min*(1-PAD); DO_hi = r_max*(1+PAD);
    DOstar = 0.5*(DO_lo+DO_hi);           % 也可改成靠近 DO_lo
else
    DOstar = DO_MANUAL; DO_lo = DO_MANUAL; DO_hi = DO_MANUAL;
end

fprintf('kz=%.2f：DO window (continuous) = [%.6f, %.6f]; choose DO*=%.6f\n', KZ_PLANE, r_min, r_max, DOstar);

%% ===== 解節點：沿 M=0 等值線解 DO^2 = r(k)^2 =====
f = @(kx) DOstar^2 - R_of_kx(kx).^2;
roots_Q1 = [];
for s = 1:size(segs,1)
    a = segs(s,1); b = segs(s,2);
    ks_seg = linspace(a, b, max(1000,round(N_SCAN/numel(segs)))+1);
    g = f(ks_seg);
    ok = ~isnan(g);
    ks_seg = ks_seg(ok); g = g(ok);
    chg = find(g(1:end-1).*g(2:end) <= 0);
    for m = 1:numel(chg)
        left = ks_seg(chg(m)); right = ks_seg(chg(m)+1);
        try
            kx_root = fzero(f, [left, right]);
            if kx_root>=a-1e-12 && kx_root<=b+1e-12
                roots_Q1(end+1,1) = kx_root; %#ok<AGROW>
            end
        catch, end
    end
end
roots_Q1 = unique(round(roots_Q1*1e9)/1e9);
nodes_Q1 = [];
for j = 1:numel(roots_Q1)
    kx = roots_Q1(j); ky = ky_from_kx(kx);
    nodes_Q1(end+1,:) = [kx, ky, KZ_PLANE]; %#ok<AGROW>
end
% 四象限複製（kz 固定）
nodes = [];
R4 = [1 1; -1 1; 1 -1; -1 -1];
for r=1:4
    nodes = [nodes; R4(r,1)*nodes_Q1(:,1), R4(r,2)*nodes_Q1(:,2), KZ_PLANE*ones(size(nodes_Q1,1),1)]; %#ok<AGROW>
end
nodes = unique(round(nodes*1e6)/1e6,'rows','stable');

if isempty(nodes)
    fprintf('找不到節點；DO* 可能不在 [%.6f, %.6f]，或需加密 N_SCAN。\n', r_min, r_max);
else
    fprintf('節點（kz=%.2f）：\n', KZ_PLANE); disp(nodes);
end

%% ===== 小盒子 3D 驗證（封閉公式，無 eig） =====
if ~isempty(nodes)
    k0 = nodes(1,:);
    [gmin3, kmin3] = local3D_verify_ftn58(k0, DOstar, A1,A2x,A2y,B1,B2x,B2y,M);
    fprintf('3D 驗證：min gap ≈ %.3e at k=(%.5f, %.5f, %.5f)\n', gmin3, kmin3);
end

%% ===== 複製到 ftn58 的參數 =====
A2_eff = 0.5*(A2x+A2y);  B2_eff = 0.5*(B2x+B2y);
fprintf('\n=== 貼回 ftn58 ===\n');
fprintf('A1=%.4f; A2=%.4f;  %% 你的腳本若只支援單一 A2\n', A1, A2_eff);
fprintf('B1=%.4f; B2=%.4f;  %% 或拆 B2x=%.4f; B2y=%.4f;\n', B1, B2_eff, B2x, B2y);
fprintf('M=%.4f; DO=%.6f; Delta=0; R=0;  %% 同步更新 DO\n', M, DOstar);

% 等向警示
if abs(B2x-B2y)<1e-12 && abs(A2x-A2y)<1e-12
    DOc = A2x*sqrt(-M/B2x); % 等向時只在單一 DO 出 nodal ring
    fprintf('\n[提醒] 等向 (B2x=B2y, A2x=A2y) 時，只有 DO=%.6f 出現 nodal ring，無離散 Weyl。\n', DOc);
end
end % main


%% ===== helpers =====
function x = clamp(x,a,b), x = min(max(x,a),b); end
function segs = contiguous_segments(xs)
if isempty(xs), segs=[]; return; end
dx = diff(xs); step = median(dx); if step==0, step = max(dx); end
brk = find(dx > 5*step);
idx = [0, brk, numel(xs)-1];
segs = zeros(numel(idx)-1,2);
for i=1:numel(idx)-1
    a = xs(idx(i)+1); b = xs(idx(i+1)+1);
    if b>a, segs(i,:) = [a,b]; end
end
segs = segs(any(segs,2),:);
end

function [gmin, kmin] = local3D_verify_ftn58(k0, DO, A1,A2x,A2y,B1,B2x,B2y,M)
% gap = 2*min_s sqrt( M_lat^2 + (A1 sin kz)^2 + (DO + s*r)^2 )
kx_loc = linspace(k0(1)-0.06, k0(1)+0.06, 25);
ky_loc = linspace(k0(2)-0.06, k0(2)+0.06, 25);
kz_loc = linspace(k0(3)-0.15*pi, k0(3)+0.15*pi, 25);
[KX,KY,KZ] = ndgrid(kx_loc,ky_loc,kz_loc);

SX = sin(KX); CX = cos(KX);
SY = sin(KY); CY = cos(KY);
SZ = sin(KZ); CZ = cos(KZ);

Mlat = M + 2*B1*(CZ-1) + 2*B2x*(CX-1) + 2*B2y*(CY-1);
R    = hypot(A2x*SX, A2y*SY);

g1 = sqrt( Mlat.^2 + (A1*SZ).^2 + (DO + R).^2 );
g2 = sqrt( Mlat.^2 + (A1*SZ).^2 + (DO - R).^2 );
gap = 2*min(g1,g2);

[gmin, idx] = min(gap(:));
[i,j,k] = ind2sub(size(gap), idx);
kmin = [KX(i,j,k), KY(i,j,k), KZ(i,j,k)];
end
