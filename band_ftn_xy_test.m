%%%          Band Structure Plot for ftn58sparse         %%%
%%% The spin polariztion is implemented in this function %%%
%%% 3/6/2016 Hans                                        %%%
%%% ---------------------------------------------------- %%%
clear all

use_tilt = 0;
T = [0.0, 0.5, 0.0]; % <---- 這邊調整傾斜

E_range = [-2 2];
Ef      = 0;%-16;

plot2DBand = 0;


kz0     = 0.00/(2*pi);   % <---- 你要固定的 kz

% 參數：費米窗大小（等能曲線粗細敏感）
% 現在 sigma 可以是多個值，程式會依序繪圖
sigmas = [0.005, 0.002, 0.0005];  % eV，可自行調整多個 sigma 值

%% Initial info. %%
tic
load ftn58sparse.mat
fprintf('載入 ftn58sparse.mat 耗時: %.2f 秒\n', toc);
norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
dd   = ftn58sparse.dd;
tt   = ftn58sparse.tt;
BR   = ftn58sparse.BR;
Sz   = [1 0;0 -1];

%% Kpoints %%%

% 2D 模式的 kpoints 網格設定
fprintf('進入 2D 模式，開始生成 kpoints 網格...\n');
N1 = 1001; N2 = 1001;                 % 解析度
kx_span = linspace(-0.5, 0.5, N1);% 掃描範圍（約化座標）
ky_span = linspace(-0.5, 0.5, N2);% 掃描範圍（約化座標）
[KX, KY] = ndgrid(kx_span, ky_span);

tic % 計時開始：生成 kpoints
% 使用向量化方式生成 kpoints
kpoints = [KX(:), KY(:), repmat(kz0, N1*N2, 1)] * 2*pi;
fprintf('向量化生成 kpoints 耗時: %.2f 秒\n', toc);

nks = size(kpoints, 1); % 更新 nks
fprintf('總共有 %d 個 k 點要計算。\n', nks);


%% Eigenvalue and Eigenvector calculations %%%
tic
% eigvec = cell(size(kpoints,1),1); % 在此模式下不再需要儲存 eigvec
nks_cal    = length(kpoints);
Ek     = zeros(nks_cal,norb);
% SPz    = zeros(nks_cal,norb); % 在此模式下不再需要計算 SPz
for ik=1:nks_cal
    time_init=tic;
    Hsparse = sparse(ii,jj,exp(1i*dd*kpoints(ik,:)').*tt,norb,norb);
    HH      = full(Hsparse);
    HH      = (HH+HH')/2 ;
    [~, Etemp] = eig(HH); % 不再需要 vec

    Ek(ik,:)  = diag(Etemp);
        %% TILT: 對每個 k 將能譜整體平移 (T·k)
    if use_tilt
        tilt_k = real(kpoints(ik,:)*T(:));   % 標準點積
        Ek(ik,:) = Ek(ik,:) + tilt_k;
    end
    %% -----
    % SSPz      = vec'*kron(Sz,eye(norb/2))*vec; % 1D 功能已移除，此段不再需要
    % SSPz      = (SSPz+SSPz')/2;
    % SPz(ik,:) = diag(SSPz);

    if mod(ik,1e4)==0
        fprintf('%3i/%i: %.3fs [%.2fm]\n',ik,nks_cal,toc,toc(time_init)/60);
    end

end
toc

% 依 norb 自動決定 reshape 的第三維
Ekk = reshape(Ek,[N1, N2, norb]);

% 使用所有能帶來找最接近 Ef 的能量距離 (這部分在 sigma 迴圈外計算，以節省重複計算)
% Ekk: [N1,N2,norb], Ef 已有；kx_span, ky_span, kz0 已有
AbsE = abs(Ekk - Ef);                     % 距離 Ef
[minDistEachK, whichBand] = min(AbsE, [], 3);
[minDist, lin] = min(minDistEachK(:));    % 全域最小
[i, j] = ind2sub([N1, N2], lin);

k_min_frac = [kx_span(i), ky_span(j), kz0];
k_min_cart = k_min_frac * 2*pi;

fprintf('[min |E-Ef|] %.6g eV @ frac (%.6f, %.6f, %.6f), band %d\n', ...
    minDist, k_min_frac, whichBand(i,j));              % [N1,N2,nb]
minAbs = min(AbsE, [], 3);              % 最接近 Ef 的距離
kz0r = kz0 *2*pi;


% --- 開始對不同的 sigma 值進行繪圖 ---
for current_sigma = sigmas
    if plot2DBand == 1
    % ========= (A) 2D-Band: kx–ky 能帶曲面（可視需要開關） =========
    % 每次迴圈都創建一個新的 figure 視窗
    figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
           'papertype','usletter','numbertitle','off','name',sprintf('kx_ky_sigma_%.4f',current_sigma), ...
           'PaperPositionMode','manual','paperorientation','portrait', ...
           'color','w');

    bands_to_plot = [2, 3];  % 例：靠近費米能的兩條
    for iband = bands_to_plot
        Z = Ekk(:,:,iband) - Ef;
        A = (abs(Z) - min(abs(Z(:)))) ./ (max(abs(Z(:))) - min(abs(Z(:))) + eps);
        surf(KX, KY, Z, ...
             'FaceColor','interp', ...
             'FaceAlpha','interp', ...
             'AlphaData', A, ...
             'AlphaDataMapping','none', ...
             'EdgeColor','none');
        hold on
    end
    lighting gouraud; box on;
    ax = gca; ax.XGrid='on'; ax.YGrid='on'; ax.ZGrid='on';
    ax.DataAspectRatio = [1 1 15];
    axis([kx_span(1) kx_span(end) ky_span(1) ky_span(end) E_range(1) E_range(2)]);
    set(gca,'FontSize',16,'TickLabelInterpreter','latex','FontWeight','bold');
    xlabel('\bf{k$_{x}$($2\pi/a$)}','interpreter','LaTex','FontSize',20,'FontWeight','bold');
    ylabel('\bf{k$_{y}$($2\pi/a$)}','interpreter','LaTex','FontSize',20,'FontWeight','bold');
    zlabel('\bf{Energy (eV)}','interpreter','LaTex','FontSize',20,'FontWeight','bold');
    colormap(parula);
    title(sprintf('kx-ky Band Surface (sigma=%.4f)', current_sigma), ...
          'Interpreter','latex','FontSize',18,'FontWeight','bold');
    end
    % ========= (B) Fermi 等能輪廓：|E - Ef| < sigma =========
    % 每次迴圈都創建一個新的 figure 視窗
    figure('position',[1020 0 760 640],'color','w','name',sprintf('fermi_contour_kxky_sigma_%.4f',current_sigma));
    contour(kx_span, ky_span, minAbs.', [current_sigma current_sigma], 'k-', 'LineWidth', 1.5);
    axis image; set(gca,'YDir','normal');
    xlabel('$k_x\,(2\pi/a)$','Interpreter','latex','FontSize',18);
    ylabel('$k_y\,(2\pi/a)$','Interpreter','latex','FontSize',18);
    title(sprintf('Fermi iso-curve: $|E-E_F|<%.5f\\,\\mathrm{eV}$, $k_z=%.3f$', ...
                  current_sigma, kz0r), ...
          'Interpreter','latex','FontSize',18,'FontWeight','bold');
end % sigma 迴圈結束