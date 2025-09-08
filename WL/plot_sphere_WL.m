% plot_S2_WL.m
clear; clc;
load BP_S2.mat  % 需要 Chern_S2, phases(Ntheta x nVB), theta, phi

% x 軸類似你原本的 axisX：把 theta 正規化到 [0,1]
x = (theta - theta(1)) / (theta(end) - theta(1));

% ------ (1) WL winding：兩條 eigenphase 線，藍點，y ∈ [-1,1] (phase/π) ------
% 把相位從 [0, 2π) 映到 [-π, π]
ph = phases;                          % Ntheta x nVB
ph(ph > pi) = ph(ph > pi) - 2*pi;     % -> [-π, π]
y1 = ph(:,1)/pi;  % 轉到 [-1,1]
y2 = ph(:,2)/pi;

figure;
hold on;
plot(x, y1, 'b.', 'MarkerSize', 6);
plot(x, y2, 'b.', 'MarkerSize', 6);
plot([0,1],[0,0],'k--');
xlabel('\theta (normalized)', 'FontSize', 14);
ylabel('phase / \pi', 'FontSize', 14);
title('WL on S^2 ( \phi-loop eigenphases vs \theta )', 'FontSize', 14);
xlim([0 1]);
ylim([-1 1]);
grid on;

% ------ (2) 累積 Chern：用相位總和的 unwrap，沿 θ 累積 ------
% 對每個本徵相位分別做 unwrap，再求和（或直接 sum 再 unwrap；兩者最後應一致）
ph_unw_1 = unwrap(phases(:,1));   % radians, 隨 θ 平滑展開
ph_unw_2 = unwrap(phases(:,2));
sum_unw  = ph_unw_1 + ph_unw_2;   % U(1) Wilson loop 的總相位（det 的相位）

% 以 θ 的起點為 0，做「累積 Chern」
cum_chern = (sum_unw - sum_unw(1)) / (2*pi);   % 應該在最後接近整數 Chern_S2

figure;
hold on;
plot(x, cum_chern, 'b.', 'MarkerSize', 6);
plot([0,1],[0,0],'k--');
xlabel('\theta (normalized)', 'FontSize', 14);
ylabel('Chern (cumulative)', 'FontSize', 14);
title(sprintf('Chern accumulation on S^2   (final \\approx %d)', Chern_S2), 'FontSize', 14);
xlim([0 1]);
% y 軸自動即可；想定 [-1,1] 也可：ylim([-1 1]);
grid on;

% 小提示：檢查終點與整數的差
final_est = cum_chern(end);
fprintf('Chern_S2 (from FHS) = %d\n', Chern_S2);
fprintf('Chern from WL winding (end) ~ %.6f (diff = %.2e)\n', final_est, final_est - Chern_S2);
