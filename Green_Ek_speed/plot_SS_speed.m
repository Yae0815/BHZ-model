clear all;
load ek.mat
%load band.mat
E_range = [-0.5 0.5];
 wef=w-EF;
 %%
 %At=Ab*5+As;
 kd=linspace(0,nk,nk);
 figure; % bulk+surface
 hold on
 %pcolor(kd,wef,(As+Ab)/pi),shading interp
 pcolor(kd,wef,(As)/pi),shading interp
% pcolor(-kd,w,At/pi),shading interp
%colorbar;
%caxis([-200 200]);
%colormap(blue_red);
% ===== 在你的圖上疊群速度箭頭：v ∝ dE/dk =====
[NE, NK] = size(As);                 % As: NE x NK（能量 × k）
[~, idx] = max(As, [], 1);           % 各 k 欄位的譜峰位置（索引）
Epk = wef(idx);                      % 對應的峰值能量 E*(k) 初估 (1 x NK)

% --- 三點拋物線微調峰值位置（假設 wef 等距）---
if NE >= 3
    dE = wef(2) - wef(1);
    for j = 2:NK-1
        i = idx(j);
        if i > 1 && i < NE
            y1 = As(i-1,j); y2 = As(i,j); y3 = As(i+1,j);
            % 頂點偏移（單位：格點）；避免除 0
            denom = (y1 - 2*y2 + y3);
            if abs(denom) > 1e-14
                delta = 0.5*(y1 - y3)/denom;  % -1~+1 之間常見
                Epk(j) = wef(i) + delta*dE;   % 微調後的 E*(k)
            end
        end
    end
end

% --- 估群速度方向：dE/dk（用 kd 的實際座標）---
dEdk = gradient(Epk, kd);            % 同單位：eV / (你的 k 軸單位)

% --- 只在 Fermi-arc 亮區畫箭頭，避免 bulk 背景 ---
linind = sub2ind([NE, NK], idx(:), (1:NK).');  % 每列的最大值位置
Amax   = As(linind);                            % 弧線的亮度估計
thr    = prctile(Amax, 98);                     % 亮度門檻（可依畫面調整 95~99）
mask   = Amax > thr;

% --- 稀疏取樣 + 箭頭長度設定 ---
step = max(4, round(NK/30));        % 箭頭密度
sel  = find(mask).';                 % 只取亮區
sel  = sel(1:step:numel(sel));       % 稀疏化
if isempty(sel), sel = 1:step:NK; end

dk   = kd(2) - kd(1);
L    = 0.6 * dk;                     % 箭頭總長（k 軸單位，可調）

vx = ones(size(sel));                % 切向向量 (Δk, ΔE) ∝ (1, dE/dk)
vy = dEdk(sel);
normv = sqrt(vx.^2 + vy.^2) + 1e-12; % 正規化避免除 0
U = L * vx ./ normv;                 % k 方向分量（資料座標）
V = L * vy ./ normv;                 % E 方向分量（資料座標）

% --- 疊加箭頭（scale=0 => U,V 已是絕對長度）---
quiver(kd(sel), Epk(sel), U, V, 0, 'w', 'LineWidth', 1.2, ...
       'MaxHeadSize', 1.5);

%（可選）標記兩個 Weyl 投影的大致位置
% xline(kd(iL), '--w', 'Alpha', 0.25);
% xline(kd(iR), '--w', 'Alpha', 0.25);

% ---- color map & scale ---------------------------------
colormap(viridis)          
caxis([0 10]);             % 0 – 50 (1/eV) 原設定

cb = colorbar;             % 新增顏色條
cb.Label.String = 'Spectral weight $A_s/\pi$ [1/eV]';
cb.Label.Interpreter = 'latex';
cb.TickDirection = 'out';
cb.FontSize = 14;
% ---------------------------------------------------------

xlim([1 1000]);
%ylim([-0.5 0.5]);
ylim(E_range);
list    = [1 501 1000];
%plot(1:nks,Ek-EF,'r-','LineWidth',1);
ylabel('\bf{Energy (eV)}','FontSize',20,'interpreter','LaTex');
ax = gca;
ax.TickDir    = 'in';
ax.FontSize   = 18;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.XTick      = list(:);
ax.YTick      = E_range(1):0.5:E_range(2);

%xticklabels({'$\mathbf{-\bar{Z}}$', '$\mathbf{\bar{\Gamma}}$', '$\mathbf{\bar{Z}}$'})
xticklabels({' ', ' ', ' '})
% ax.XTickLabel = {};
ax.TickLabelInterpreter='latex';