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