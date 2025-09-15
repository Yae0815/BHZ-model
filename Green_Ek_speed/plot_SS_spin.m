clear all;
load ek.mat
%load band.mat
comp = 'y';
E_range = [-0.3 0.3];
 wef=w-EF;
 %%
 %At=Ab*5+As;
 kd=linspace(0,nk,nk);
 figure; % bulk+surface
 hold on
 if comp == 'x'
 pcolor(kd,wef,Asx/pi),shading interp
 end
 if comp == 'y'
 pcolor(kd,wef,Asy/pi),shading interp
 end
 if comp == 'z'
 pcolor(kd,wef,Asz/pi),shading interp
 end
% pcolor(-kd,w,At/pi),shading interp
colorbar;
caxis([-5 5]);
colormap(blue_red);
cb = colorbar;
%caxis([0 2]);
xlim([0 1000]);
ylim([-0.3 0.3]);
list    = [1 501 1000];
% plot(1:nks,Ek-EF,'k-','LineWidth',1);
ylabel('\bf{Energy (eV)}','FontSize',24,'interpreter','LaTex');
cb.Label.String = sprintf('Spin Spectral Function $A_{s,%c}/\\pi$ [1/eV]', comp);
cb.Label.Interpreter = 'latex';
ax = gca;
ax.TickDir    = 'in';
ax.FontSize   = 18;
ax.FontWeight = 'bold';
ax.TickLength = [0.02 0.02];
ax.XTick      = list(:);
ax.YTick      = E_range(1):1:E_range(2);
%xticklabels({'$\mathbf{-\bar{Z}}$', '$\mathbf{\bar{\Gamma}}$', '$\mathbf{\bar{Z}}$'})
%xticklabels({'$-0.1 \pi$', '$\mathbf{\bar{\Gamma}}$', '$0.1 \pi$'})
xticklabels({' ', ' ', ' '})
% ax.XTickLabel = {'\bf{Y}' '\bf{$\Gamma$}' '\bf{-Y}'};
% ax.XTickLabel = {};
ax.TickLabelInterpreter='latex';