%%%          Band Structure Plot for ftn58sparse         %%%
%%% The spin polariztion is implemented in this function %%%
%%% 3/6/2016 Hans                                        %%%
%%% ---------------------------------------------------- %%%

E_range = [-1 1];
Ef      = 0;%-16;
isSP    = 1;

%% Initial info. %%
load Sftn58sparse.mat
ftn58sparse = Sftn58sparse;
norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
dd   = ftn58sparse.dd;
tt   = ftn58sparse.tt;
BR   = ftn58sparse.BR2D;
Sx = [0 1; 1 0];
Sy = [0 -1i; 1i 0];
Sz = [1 0; 0 -1];

%% Kpoints %%%
nk = 100;

n = 0.5;
kpath = 'x';  % <----改這裡選路徑

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


list    = [1 101 200];
kpoints = [p1(1:nk,:);p2(1:nk,:)]*2*pi;



function h = cplot(x, y, c, varargin)
% CPLOT  Color-graded 2-D line: x, y, scalar c
x = x(:); y = y(:); c = c(:);          % 確保列向量
x = [x; NaN]; y = [y; NaN]; c = [c; NaN];  % NaN 斷線
h = patch('XData',x,'YData',y, ...
          'CData',c,'EdgeColor','interp', ...
          'FaceColor','none', varargin{:});
end



    

%% Eigenvalue and Eigenvector calculations %%%
tic
eigvec = cell(size(kpoints,1),1);
nks    = length(kpoints);
Ek     = zeros(nks,norb);
SPx = zeros(nks,norb);   %〈Sx〉
SPy = zeros(nks,norb);   %〈Sy〉
SPz = zeros(nks,norb);   %〈Sz〉

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
    SSPx = vec' * kron(Sx, eye(norb/2)) * vec;  SSPx = (SSPx+SSPx')/2;
    SSPy = vec' * kron(Sy, eye(norb/2)) * vec;  SSPy = (SSPy+SSPy')/2;
    SSPz = vec' * kron(Sz, eye(norb/2)) * vec;  SSPz = (SSPz+SSPz')/2;

    SPx(ik,:) = diag(SSPx);
    SPy(ik,:) = diag(SSPy);
    SPz(ik,:) = diag(SSPz);
    if mod(ik,1e3)==0
        fprintf('%3i/%i: %.3fs [%.2fm]\n',ik,nks,toc,toc(time_init)/60);
    end
    
end
toc


[~, Emin] = min(abs(Ek(1,:)-Ef));
disp(Emin);

function gap = direct_gap(Ek_row, Ef)
    Evbm = max(Ek_row(Ek_row <= Ef));   % valence-band max
    Ecbm = min(Ek_row(Ek_row >= Ef));   % conduction-band min
    gap  = Ecbm - Evbm;                % direct gap (可為負—代表交疊)
end



%% ---------- Spin-resolved Plotting ----------
if isSP
    % 產生平滑色帶（藍←0→紅），並在每張圖上套用
    try
        cm = blue_red(256);     % 自訂函數：請確保 blue_red.m 在 path 中
    catch
        error('找不到 blue_red.m，請將它放到 MATLAB path。');
    end
    spmap = flipud(cm);         % 翻轉：低值→藍，高值→紅

    SPcell   = {SPx, SPy, SPz};
    slabel   = {'S_x','S_y','S_z'};

    for ic = 1:3
        figure('Color','w');  hold on

        %─── 能帶着色 ───
        for ii = 1:norb
            cplot(1:nks, Ek(:,ii)-Ef, SPcell{ic}(:,ii), 'LineWidth',1);
        end

        %─── 每張圖套用 colormap 與色階 ───
        colormap(spmap);
        caxis([-1 1]);          % 若要自動縮放，可改成 caxis(max(abs(SPcell{ic}(:)))*[-1 1])

        %─── 色條 ───
        cb                      = colorbar;
        cb.TickLabelInterpreter = 'latex';
        cb.FontSize             = 18;
        cb.Label.Interpreter    = 'latex';
        cb.Label.FontSize       = 20;
        cb.Label.String         = ['$\langle ', slabel{ic}, ' \rangle$'];

        %─── 軸外觀與標題 ───
        decorate_axis(list, E_range, nks, kpath);
        title(['Spin-Resolved Bands : $', slabel{ic}, '$'], ...
              'Interpreter','latex', 'FontSize',18, 'FontWeight','bold');
    end
else
    ...
end



%% ---------- 軸設定函式 ----------

function decorate_axis(list, E_range, nks, kpath)
    % 垂直線：Γ 及 ±Y
    for il = 2:numel(list)-1
        line([list(il) list(il)], E_range, ...
             'LineStyle','-', 'LineWidth',0.1, 'Color','k');
    end
    % 0 eV 基準線
    line([0 nks], [0 0], ...
         'LineStyle','--', 'LineWidth',0.5, 'Color','k');

    % 軸格式
    axis([1 nks E_range]);
    ax                       = gca;
    ax.TickDir               = 'in';
    ax.FontSize              = 14;
    ax.FontWeight            = 'bold';
    ax.TickLength            = [0.02 0.02];
    ax.XTick                 = list(:);
    ax.YTick                 = E_range(1):1:E_range(2);    
    if kpath == 'x'
        ax.XTickLabel            = {'\bf{-X}','\bf{$\Gamma$}','\bf{X}'};
    end    
    if kpath == 'y'
        ax.XTickLabel            = {'\bf{-Y}','\bf{$\Gamma$}','\bf{Y}'};
    end
    if kpath == 'z'
        ax.XTickLabel            = {'\bf{-Z}','\bf{$\Gamma$}','\bf{Z}'};
    end
    ax.LineWidth             = 1.2;
    ax.TickLabelInterpreter  = 'latex';
    ylabel('\bf{Energy (eV)}', ...
           'Interpreter','latex', 'FontSize',24, 'FontWeight','bold');
end