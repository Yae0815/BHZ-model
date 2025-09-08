% ==============================================================
% batch_surface_spin.m
% 跑 top/bot × x/y/z，並自動存成 top_kx_sx.jpg 等
% ==============================================================

clear; clc;

% 你的原繪圖腳本檔名（請改成實際檔名）
target_script = 'half_slab_spin.m';

% 輸出資料夾
outdir = 'figs';
if ~exist(outdir,'dir'); mkdir(outdir); end

% 要掃的參數
sides  = {'top','bot'};
kpaths = {'x','y','z'};

% （可選）你也能在這邊統一覆蓋原腳本的參數，例如：
% use_surface_overlap = 1;
% layers_take = 2;
% norb_per_layer_nospin = 4;
% M_keep = 8;
% E_range = [-1 1];
% Ef = 0;
% isSP = 1;
% nk = 100;

for iside = 1:numel(sides)
    for ikp = 1:numel(kpaths)
        % 把變數丟到 base workspace，讓 script 接到
        side      = sides{iside};
        kpath     = kpaths{ikp};
        auto_save = true;                % 觸發自動存圖
        % outdir 已在本檔定義，同名變數會被原檔讀到
        fprintf('>>> Running %s | side=%s, kpath=%s\n', target_script, side, kpath);

        % 直接執行原腳本（原腳本開頭有 exist 檢查就會採用這裡的設定）
        run(target_script);
    end
end

fprintf('All done. Files saved under: %s\n', outdir);
