%% === Chern numbers at kz = 0 and kz = pi from Wilson-loop phases ===
clear; load BP.mat   % 需含: phase [nVB, nky, nkz], ind, axisX (可選)
% 說明：phase 存的是「以 π 為單位」的相位；下方會轉回弧度計算

[nVB, nky, nkz] = size(phase);
assert(nkz>=2, 'BP.mat 需至少含 kz=0 與 kz=pi 兩個切片');

chern = zeros(1,nkz);
craw  = zeros(1,nkz);
err   = zeros(1,nkz);

for ikz = 1:nkz
    % 取出該 kz 的相位並轉回弧度
    ph = squeeze(phase(:,:,ikz))*pi;   % [nVB, nky], radians

    % —— 連續帶配對（沿 ky）——
    ph_sorted = zeros(size(ph));
    ph_sorted(:,1) = sort(ph(:,1));    % 起點先排序

    for iy = 2:nky
        prev = ph_sorted(:,iy-1);
        curr = ph(:,iy);

        % 最短圓距離成本矩陣（考慮 2π 等價）
        cost = zeros(nVB);
        for a = 1:nVB
            for b = 1:nVB
                d = curr(b) - prev(a);
                d = atan2(sin(d), cos(d));   % 映到 (-pi,pi]
                cost(a,b) = abs(d);
            end
        end

        % 簡單貪婪指派（nVB 小時穩定）
        used = false(1,nVB);
        curr_sorted = zeros(nVB,1);
        for a = 1:nVB
            [~, idx] = min(cost(a,:) + 1e3*used);
            used(idx) = true;
            curr_sorted(a) = curr(idx);
        end
        ph_sorted(:,iy) = curr_sorted;
    end

    % —— 展相並計算總繞行數（= Chern）——
    ph_unw = unwrap(ph_sorted, [], 2);      % 對 ky 展相
    dphi   = ph_unw(:,end) - ph_unw(:,1);   % 每條佔據帶淨變化
    C_raw  = sum(dphi)/(2*pi);              % 合併後除以 2π
    Chern  = round(C_raw);                  % 取最近整數
    chern(ikz) = Chern; craw(ikz) = C_raw; err(ikz) = C_raw - Chern;

    fprintf('kz index %d  ->  Chern = %d   (raw = %.6f, deviation = %.2e)\n', ...
            ikz, Chern, C_raw, err(ikz));

    % 視覺檢查：畫出展相後的相位流（以 π 為單位顯示）
    figure; hold on;
    plot(1:nky, (ph_unw'/pi), 'b.');
    yline(0,'k--');
    xlabel('ky index'); ylabel('phase / \pi');
    if isfield(ind,'kz')
        title(sprintf('Unwrapped Wilson phases @ kz=%.3f | Chern=%d', ind.kz(ikz), Chern));
    else
        title(sprintf('Unwrapped Wilson phases @ kz index %d | Chern=%d', ikz, Chern));
    end
end

% 便捷讀數
if nkz >= 2
    fprintf('\nSummary:\n  Chern(kz=0) = %d\n  Chern(kz=pi) = %d\n', chern(1), chern(2));
end

