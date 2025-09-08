%%%          Band Structure Plot for ftn58sparse         %%%
%%% The spin polariztion is implemented in this function %%%
%%% 3/6/2016 Hans                                        %%%
%%% ---------------------------------------------------- %%%

E_range = [-1 1];
Ef      = 0;%-16;
isSP    = 1;
is2D    = 0;

%% Initial info. %%
load Sftn58sparse.mat
ftn58sparse = Sftn58sparse;
norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
dd   = ftn58sparse.dd;
tt   = ftn58sparse.tt;
BR   = ftn58sparse.BR2D;
Sz   = [1 0;0 -1];

%% Kpoints %%%
nk = 100;

n = 0.5;
p1 = [linspace(0.0,0.0,nk+1)' linspace(0.5,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
p2 = [linspace(0.0,0.0,nk+1)' linspace(0.0,-0.5,nk+1)' linspace(0.0,0.0,nk+1)'];
%p1 = [linspace(0.5,0.0,nk+1)' linspace(0.0,0.0,nk+1)' linspace(0.0,0.0,nk+1)'];
%p2 = [linspace(0.0,-0.5,nk+1)' linspace(0.0,0.0,nk+1)' linspace(0.0,0.0,nk+1)']

list    = [1 101 200];
kpoints = [p1(1:nk,:);p2(1:nk,:)]*2*pi;

% 生成 path 後馬上記下索引
idxY      = 1;          % k = (0, 0.5, 0)
idxGamma  = nk+1;       % k = (0, 0,   0)  ── list(2)
idxGamma_near3  = nk+4;
idxGamma_near6  = nk+7;
idxmY     = 2*nk;     % k = (0,-0.5, 0)


function h = cplot(x, y, c, varargin)
% CPLOT  Color-graded 2-D line: x, y, scalar c
x = x(:); y = y(:); c = c(:);          % 確保列向量
x = [x; NaN]; y = [y; NaN]; c = [c; NaN];  % NaN 斷線
h = patch('XData',x,'YData',y, ...
          'CData',c,'EdgeColor','interp', ...
          'FaceColor','none', varargin{:});
end


if is2D==1
    N1 = 251; N2 = 251;
    p1 = linspace(-0.0,0.2,N1); 
    p2 = linspace(-0.0,0.2,N2);  
    [K1,K2] = meshgrid(p1,p2);

    kk = 1;
    for n1=1:N1
        for n2=1:N2
            kpoints(kk,:) = [p1(n1) p2(n2) 0]*2*pi;
            kk = kk + 1;
        end
    end

    nks     = length(kpoints);
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
    SSPz      = vec'*kron(Sz,eye(norb/2))*vec;
    SSPz      = (SSPz+SSPz')/2;
    SPz(ik,:) = diag(SSPz); 
    
    if mod(ik,1e3)==0
        fprintf('%3i/%i: %.3fs [%.2fm]\n',ik,nks,toc,toc(time_init)/60);
    end
    
end
toc


if is2D==1
    Ekk = reshape(Ek,[N1,N2,4]);

    %% 2D-Band
    figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
           'papertype','usletter','numbertitle','off','name','kx_kz',...
           'PaperPositionMode','manual','paperorientation','portrait',...
           'color','w');

    for iband = [2 3]
    %     surface(K1,K2,Ekk(:,:,iband)-Ef,'LineStyle','-','LineWidth',0.01,...
    %             'EdgeColor','none');
        surf(K1,K2,Ekk(:,:,iband)-Ef,'FaceAlpha','interp',...
            'FaceColor','blue',...    
            'AlphaDataMapping','scaled',...
            'AlphaData',gradient(abs(Ekk(:,:,iband)-Ef))*1e2,...        
            'EdgeColor','none');
        hold on
    end

    lighting gouraud;
    box on;
    ax = gca;
    ax.XGrid= 'on';
    ax.YGrid= 'on';
    ax.ZGrid= 'on';
    ax.DataAspectRatio = [1 1 15];
    axis([K1(1) K1(end) K2(1) K2(end) -1 1]);

    set(gca,'FontSize',16,'TickLabelInterpreter','latex','FontWeight','bold');
    xlabel('\bf{k$_{x}$($2\pi/a$)}','interpreter','LaTex','FontSize',20,'FontWeight','bold');
    ylabel('\bf{k$_{y}$($2\pi/a$)}','interpreter','LaTex','FontSize',20,'FontWeight','bold');
    zlabel('\bf{Energy (eV)}','interpreter','LaTex','FontSize',20,'FontWeight','bold');
    return
end

[~, Emin] = min(abs(Ek(1,:)-Ef));
disp(Emin);

function gap = direct_gap(Ek_row, Ef)
    Evbm = max(Ek_row(Ek_row <= Ef));   % valence-band max
    Ecbm = min(Ek_row(Ek_row >= Ef));   % conduction-band min
    gap  = Ecbm - Evbm;                % direct gap (可為負—代表交疊)
end

gap_G  = direct_gap(Ek(idxGamma,:), Ef);
gap_G3  = direct_gap(Ek(idxGamma_near3,:), Ef);
gap_G6 = direct_gap(Ek(idxGamma_near6,:), Ef);
gap_G  = direct_gap(Ek(idxGamma,:), Ef);
gap_Y  = direct_gap(Ek(idxY     ,:), Ef);
gap_mY = direct_gap(Ek(idxmY   ,:), Ef);
fprintf('\nDirect gaps (eV): Γ = %.4f | Γ+3 = %.4f | Γ+6 = %.4f\n',...
         gap_G, gap_G3, gap_G6);
fprintf('\nDirect gaps(Near by Gamma) (eV): Γ = %.4f | Y = %.4f | –Y = %.4f\n',...
         gap_G, gap_Y, gap_mY);


%% Plotting %%%
figure

if isSP==1
    ncolor   = 1e4;
    spmap_p  = [linspace(1,1,ncolor+1)' linspace(0,1,ncolor+1)' linspace(0,1,ncolor+1)'];
    spmap_n  = [linspace(1,0,ncolor)' linspace(1,0,ncolor)' linspace(1,1,ncolor)'];
    spmap    = [spmap_p(1:ncolor,:);spmap_n];

    for ii=1:norb
        cplot(1:nks,Ek(:,ii)-Ef,SPz(:,ii),'LineWidth',1);
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
ax.XTickLabel = {'\bf{Y}' '\bf{$\Gamma$}' '\bf{-Y}'};
% ax.XTickLabel = {'\bf{X}' '\bf{$\Gamma$}' '\bf{X}'};
% ax.XTickLabel = {};
ax.LineWidth  = 1.2;
ax.TickLabelInterpreter='latex';