function refine_from_prev_DO(DO_prev)
% Second pass: use DO_prev from your first run, solve nodes on kz plane by
% continuous contour method (no min-gap bisection), and verify in 3D.
% ftn58-consistent mass sign: Mlat = M + 2B1(cos kz - 1) + 2B2x(cos kx - 1) + 2B2y(cos ky - 1)

if nargin==0, DO_prev = 0.061012; end  % <-- 把你的第一次結果丟進來

%% ===== Parameters (match your printout) =====
A1  = 0.2;
A2x = 0.2;  A2y = 0.2;            % 你若有用 A2 異向性就改這裡
B1  = 0.3000;
B2x = 0.3000;                     % 對應 cos kx
B2y = 0.3308;                     % 對應 cos ky（你上次貼的是 0.3154 的 1.05 倍）
M   = -0.0300;                    % 注意：這裡用 ftn58 的 M
C=0; D1=0; D2=0; %#ok<NASGU>

KZ_PLANE = 0;                     % 只加 DO 的模型，節點鎖在 sin kz = 0；先檢查 kz=0
N_SCAN   = 12000;                 % 一維掃描密度（越大越準）

%% ===== 定義 kz 平面上的連續函數（ftn58 號誌）=====
cz = cos(KZ_PLANE);
Cy_from_kx = @(kx) 1 - ((-M) - 2*B1*(1-cz) - 2*B2x*(1-cos(kx)))/(2*B2y);  % M=0 ⇒ cos ky
valid_ky   = @(kx) abs(Cy_from_kx(kx))<=1+1e-14;
ky_from_kx = @(kx) acos( clamp(Cy_from_kx(kx), -1, 1) );
R_of_kx    = @(kx) hypot( A2x*sin(kx), A2y*sin(ky_from_kx(kx)) );         % r(k)=√[(A2x sin kx)^2+(A2y sin ky)^2]

%% ===== 找到 M=0 的有效區段並沿等值線解 DO_prev^2 = r(k)^2 =====
ks = linspace(0, pi, N_SCAN+1);
mask = arrayfun(valid_ky, ks);
segs = contiguous_segments(ks(mask));
if isempty(segs)
    error('kz=%.2f 平面沒有 M=0 等值線；請改 M 或試 kz=pi。', KZ_PLANE);
end

f = @(kx) DO_prev^2 - R_of_kx(kx).^2;
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
nodes = unique(round(nodes*1e6)/1e6, 'rows', 'stable');

fprintf('--- Second pass with DO = %.6f ---\n', DO_prev);
if isempty(nodes)
    fprintf('找不到節點；這通常代表 DO 還在視窗邊緣外或 M/B2 設定不吻合。\n');
else
    disp('Weyl nodes on kz=0:'); disp(nodes);
end

%% ===== 小盒子 3D 驗證（封閉公式，無 eig） =====
if ~isempty(nodes)
    k0 = nodes(1,:);  % 拿第一顆做局部檢查
    [gmin, kmin] = local3D_verify_ftn58(k0, DO_prev, A1,A2x,A2y,B1,B2x,B2y,M);
    fprintf('Local 3D check: min gap ≈ %.3e at k=(%.5f, %.5f, %.5f)\n', gmin, kmin);
end

%% ===== 供貼回 ftn58 的參數 =====
A2_eff = 0.5*(A2x+A2y);  B2_eff = 0.5*(B2x+B2y);
fprintf('\n=== 貼回 ftn58（第二回合） ===\n');
fprintf('A1=%.4f; A2=%.4f;  B1=%.4f; B2=%.4f;\n', A1, A2_eff, B1, B2_eff);
fprintf('%% 或拆：B2x=%.4f; B2y=%.4f;\n', B2x, B2y);
fprintf('M=%.4f; DO=%.6f; Delta=0; R=0;\n', M, DO_prev);

end % ===== end main =====

%% ===== helpers =====
function x = clamp(x,a,b), x = min(max(x,a),b); end
function segs = contiguous_segments(xs)
if isempty(xs), segs=[]; return; end
dx = diff(xs); step = median(dx); if step==0, step=max(dx); end
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
% gap = 2*min_s sqrt( Mlat^2 + (A1 sin kz)^2 + (DO + s*r)^2 ), with ftn58 sign
kx_loc = linspace(k0(1)-0.05, k0(1)+0.05, 25);
ky_loc = linspace(k0(2)-0.05, k0(2)+0.05, 25);
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
