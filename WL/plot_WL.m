clear all
%close all
load BP.mat
%%
% figure('position',[150 0 450 360],'paperposition',[0.25 0.25 30 10],...
%     'papertype','usletter','numbertitle','off',...
%     'PaperPositionMode','manual','paperorientation','portrait',...
%     'color','w');
%%
figure
hold on
nkx = length(ind.kx(1,:))/2+1;
plot(axisX(1,nkx:end),phase(:,nkx:end,1)','b.','MarkerSize',6)
hold on
plot([0,1],[0,0],'k--');
xlabel('ky','FontSize',14);
ylabel('phase','FontSize',14);
title('WL z=0','FontSize',14);
xlim([0 1]);
ylim([-1 1]);
%%
figure
hold on
nkx = length(ind.kx(1,:))/2+1;
plot(axisX(1,nkx:end),phase(:,nkx:end,2)','b.','MarkerSize',6)
hold on
plot([0,1],[0,0],'k--');
xlabel('ky','FontSize',14);
ylabel('phase','FontSize',14);
title('WL z=pi','FontSize',14);
xlim([0 1]);
ylim([-1 1]);
%%
figure
hold on
nkx = [1:length(ind.kx(1,:))];
plot(nkx,CkWiLp(:,1)','b.','MarkerSize',6)
hold on
plot([0,1],[0,0],'k--');
xlabel('ky','FontSize',14);
ylabel('phase','FontSize',14);
title('Chern z=0','FontSize',14);
ylim([-1 1]);
%%
figure
hold on
nkx = [1:length(ind.kx(1,:))];
plot(nkx,CkWiLp(:,2)','b.','MarkerSize',6)
hold on
plot([0,1],[0,0],'k--');
xlabel('ky','FontSize',14);
ylabel('phase','FontSize',14);
title('Chern z=pi','FontSize',14);
ylim([-1 1]);