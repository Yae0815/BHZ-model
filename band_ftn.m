
%%%          Band Structure Plot for ftn58sparse         %%%
%%% The spin polariztion is implemented in this function %%%
%%% 3/6/2016 Hans                                        %%%
%%% ---------------------------------------------------- %%%
clear all

use_tilt = 0;
T = [0.0, 0.5, 0.0]; % <---- 這邊調整傾斜

E_range = [-3 3];
Ef      = 0;%-16;
isSP    = 0;
kpath = 'x';
is2D    = 0;
kx0 = 0.0;
% --- 費米窗大小 (等能曲線) ---
sigma = 0.001;  % eV，可自行調整
%% Initial info. %%
load ftn58sparse.mat
norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
dd   = ftn58sparse.dd;
tt   = ftn58sparse.tt;
BR   = ftn58sparse.BR;
Sz   = [1 0;0 -1];

%% Kpoints %%%
nk = 400;

n = 0.5;
if kpath == 'x'
p1 = [linspace(-0.5,0.0,nk+1)' linspace(0.0,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
p2 = [linspace(0.0,0.5,nk+1)' linspace(0.0,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
end
if kpath == 'y'
p1 = [linspace(0.0,0.0,nk+1)' linspace(-0.5,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
p2 = [linspace(0.0,0.0,nk+1)' linspace(0.0,0.5,nk+1)' linspace(0.0,0.0,nk+1)'];
end
if kpath == 'z'
p1 = [linspace(0.0,0.0,nk+1)' linspace(0.0,0.0,nk+1)' linspace(-0.5,0.0,nk+1)'];
p2 = [linspace(0.0,0.0,nk+1)' linspace(0.0,0.0,nk+1)' linspace(0.0,0.5,nk+1)'];
end
if kpath == 'r'
p1 = [linspace(-0.5,0.0,nk+1)' linspace(-0.5,0.0,nk+1)' linspace(-0.5,0.0,nk+1)'];
p2 = [linspace(0.0,0.5,nk+1)' linspace(0.0,0.5,nk+1)' linspace(0.0,0.5,nk+1)'];
end
if kpath == 'm'
p1 = [linspace(0.0,0.0458,nk)' linspace(0.122,0.122,nk)' linspace(0.00,0.00,nk)'];
p2 = [linspace(0.0458,0.0916,nk)' linspace(0.122,0.122,nk)' linspace(0.00,0.00,nk)'];
end

list    = [1 401 800];
kpoints = [p1(1:nk,:);p2(1:nk,:)]*2*pi;

if is2D==1
    % 掃描 ky–kz 平面，固定 kx = kx0
    N1 = 251; N2 = 251;                 % 解析度
    ky_span = linspace(-0.2, 0.2, N1);  % 自行調整掃描範圍（約化座標）
    kz_span = linspace(-0.2, 0.2, N2);  % 自行調整掃描範圍（約化座標）
    [KY, KZ] = meshgrid(ky_span, kz_span);

    kk = 1;
    for n1=1:N1
        for n2=1:N2
            % 這裡用 [kx0, ky, kz]，並在最後乘上 2*pi（與你原程式一致）
            kpoints(kk,:) = [kx0, KY(n2,n1), KZ(n2,n1)] * 2*pi;
            kk = kk + 1;
        end
    end
    nks = length(kpoints);
end
    

%% Eigenvalue and Eigenvector calculations %%%
tic
eigvec = cell(size(kpoints,1),1);
nks    = length(kpoints);
Ek     = zeros(nks,norb);
SPz    = zeros(nks,norb);
for ik=1:nks
    time_init=tic;
    Hsparse = sparse(ii,jj,exp(1i*dd*kpoints(ik,:)').*tt,norb,norb);
    HH      = full(Hsparse);
    HH      = (HH+HH')/2 ;
    [vec, Etemp] = eig(HH);
     
    Ham{ik,1}    = HH;
    eigvec{ik,1} = vec;

%     A(ik*norb-norb+1:ik*norb,1) = ik;
%     A(ik*norb-norb+1:ik*norb,2) = diag(Etemp);
    
    Ek(ik,:)  = diag(Etemp);
        %% TILT: 對每個 k 將能譜整體平移 (T·k)
    if use_tilt
        tilt_k = real(kpoints(ik,:)*T(:));   % 標準點積
        Ek(ik,:) = Ek(ik,:) + tilt_k;
    end
    %% -----
    SSPz      = vec'*kron(Sz,eye(norb/2))*vec;
    SSPz      = (SSPz+SSPz')/2;
    SPz(ik,:) = diag(SSPz); 
    
    if mod(ik,1e3)==0
        fprintf('%3i/%i: %.3fs [%.2fm]\n',ik,nks,toc,toc(time_init)/60);
    end
    
end
toc

if is2D==1
    % 依 norb 自動決定 reshape 的第三維
    Ekk = reshape(Ek,[N1, N2, norb]);

    % ========= (A) 2D-Band: ky–kz 曲面（原功能，保留可調參） =========
    figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
           'papertype','usletter','numbertitle','off','name','ky_kz',...
           'PaperPositionMode','manual','paperorientation','portrait',...
           'color','w');

    bands_to_plot = [2, 3];  % 例：靠近費米能的兩條
    for iband = bands_to_plot
        Z = Ekk(:,:,iband) - Ef;
        A = (abs(Z) - min(abs(Z(:)))) ./ (max(abs(Z(:))) - min(abs(Z(:))) + eps);
        surf(KY, KZ, Z, ...
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
    axis([ky_span(1) ky_span(end) kz_span(1) kz_span(end) E_range(1) E_range(2)]);
    set(gca,'FontSize',16,'TickLabelInterpreter','latex','FontWeight','bold');
    xlabel('\bf{k$_{y}$($2\pi/a$)}','interpreter','LaTex','FontSize',20,'FontWeight','bold');
    ylabel('\bf{k$_{z}$($2\pi/a$)}','interpreter','LaTex','FontSize',20,'FontWeight','bold');
    zlabel('\bf{Energy (eV)}','interpreter','LaTex','FontSize',20,'FontWeight','bold');
    colormap(parula);

    % ========= (B) Fermi-slice 2D 分布：|E - Ef| < delta =========
    % 參數：費米窗與要計數的能帶範圍
    fermi_window = 0.01;              % eV，可自行調整
    bands_fermi  = [];                 % 例如指定 [Emin-1:Emin+2]；空=[] 表示使用所有能帶
    if isempty(bands_fermi)
        bands_fermi = 1:norb;
    else
        bands_fermi = bands_fermi(bands_fermi>=1 & bands_fermi<=norb);
    end

    % 計算每個 (ky,kz) 有多少能帶落在費米窗內
    AbsE = abs(Ekk(:,:,bands_fermi) - Ef);     % [N1,N2,nb]
    nHit = sum(AbsE < fermi_window, 3);        % 命中帶數（熱度）
    minAbs = min(AbsE, [], 3);                 % 最接近 Ef 的距離，用於等高線

    
    

    
    % --- 繪製等能輪廓 ---
    figure('position',[1020 0 760 640],'color','w','name','fermi_contour');
    contour(ky_span, kz_span, minAbs', [sigma sigma], 'k-', 'LineWidth', 1.5);
    
    axis image; set(gca,'YDir','normal');
    xlabel('$k_z\,(2\pi/a)$','Interpreter','latex','FontSize',18);
    ylabel('$k_y\,(2\pi/a)$','Interpreter','latex','FontSize',18);
    title(sprintf('Fermi iso-curve: $|E-E_F|<%.3f\\,\\mathrm{eV}$, $k_x=%.3f$', ...
                  sigma, kx0), ...
          'Interpreter','latex','FontSize',18,'FontWeight','bold');
    return
end

[~, Emin] = min(abs(Ek(1,:)-Ef));
disp(Emin);

%% Plotting %%%
figure

if isSP==1
    ncolor   = 1e4;
    spmap_p  = [linspace(1,1,ncolor+1)' linspace(0,1,ncolor+1)' linspace(0,1,ncolor+1)'];
    spmap_n  = [linspace(1,0,ncolor)' linspace(1,0,ncolor)' linspace(1,1,ncolor)'];
    spmap    = [spmap_p(1:ncolor,:);spmap_n];

    for ii=1:norb
        cplot(1:nks,Ek(:,ii)-Ef,SPz(:,ii),'-','LineWidth',1);
        colormap(spmap);
        hold on
    end
    colorbar('TickLabelInterpreter','LaTex','FontSize',18);
else
    for ii=1:norb
        plot(1:nks,Ek(:,ii)-Ef,'b-','LineWidth',1.5);
        hold on
    end
end

for il = 2:size(list,2)-1
line('XData', [list(il) list(il)], 'YData', [E_range(1) E_range(2)], 'LineStyle', '-', ...
    'LineWidth', 0.1, 'Color','k');
end
line('XData', [0 nks], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color','k');
ylabel('\bf{Energy (eV)}','FontSize',24,'interpreter','LaTex');

axis([1 nks E_range(1) E_range(2)]);
ax = gca;
ax.TickDir    = 'in';
ax.FontSize   = 14;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.XTick      = list(:);
ax.YTick      = E_range(1):1:E_range(2);
    if kpath == 'x'
        ax.XTickLabel = {'\bf{-X}','\bf{$\Gamma$}','\bf{X}'};
    end
    if kpath == 'y'
        ax.XTickLabel = {'\bf{-Y}','\bf{$\Gamma$}','\bf{Y}'};
    end
    if kpath == 'z'
        ax.XTickLabel = {'\bf{-Z}','\bf{$\Gamma$}','\bf{Z}'};
    end
    if kpath == 'r'
        ax.XTickLabel = {'\bf{-R}','\bf{$\Gamma$}','\bf{R}'};
    end
% ax.XTickLabel = {};
ax.LineWidth  = 1.2;
ax.TickLabelInterpreter='latex';