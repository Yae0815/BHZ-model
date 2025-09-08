function [P_top, P_bot, info] = make_surface_projectors(ftn58sparse, norb, Ns)
% 自動從 ftn58sparse 推斷 slab 層結構並建立上/下表面投影器（含自旋）
% 需要：含自旋總軌域數 norb，表面取 Ns 層
% 嘗試讀取下列任一字段作為無自旋軌域的位置 (Norb0 x 3)：
%   ftn58sparse.orbR, .Rorb, .xyz, .pos, .orb_pos（擇一存在即可）
%
% 回傳：
%   P_top, P_bot : (norb x norb) 稀疏投影器（含自旋）
%   info         : 結構，含 fields:
%                  .z (Norb0x1), .z_levels (層心 z), .layer_of_orb (Norb0x1)
%                  .NL（自動偵測的總層數）, .NorbL（每層無自旋軌域數)

    Norb0 = norb/2;   % 無自旋軌域數（假設基底為 [spin]⊗[orb(no-spin)]）
    assert(abs(Norb0-round(Norb0))<1e-12,'norb 需為偶數（含自旋）。');
    Norb0 = round(Norb0);

    % 1) 找無自旋軌域的位置陣列 (Norb0 x 3)
    cand = {'orbR','Rorb','xyz','pos','orb_pos'};
    pos = [];
    for k = 1:numel(cand)
        if isfield(ftn58sparse, cand{k})
            pos = ftn58sparse.(cand{k});
            break
        end
    end
    if isempty(pos)
        error(['找不到軌域座標。請在 ftn58sparse 中提供 (Norb0 x 3) 的位置字段（例如 orbR/xyz）。',...
               ' 目前可用鍵：orbR, Rorb, xyz, pos, orb_pos']);
    end
    if size(pos,2)~=3
        error('軌域位置陣列需為 (Norb0 x 3)。');
    end
    if size(pos,1)~=Norb0
        % 有些資料會存成「含自旋」數；若相等則自動拆半
        if size(pos,1)==norb
            % 假設「含自旋座標」是兩份相同 → 取前半
            pos = pos(1:Norb0, :);
        else
            error('位置數量與無自旋軌域數不一致，請檢查資料。');
        end
    end

    z = pos(:,3);
    % 2) 以 z 值自動分層：用自適應容忍度把相近 z 歸為同層
    z_sorted = sort(z);
    dz = diff(z_sorted);
    dz = dz(dz>eps);
    if isempty(dz)
        error('所有 z 完全相同，無法分層。請檢查座標。');
    end
    tol = 0.5*min(dz);                 % 最小間距的一半做容忍
    % 把 z 量化到格點以便 group
    zq = round(z/tol)*tol;
    [z_levels,~,layer_of_orb] = unique(zq);  % 自動給每個軌域一個 layer id
    NL = numel(z_levels);

    % 確認每層軌域數是否一致（常見 slab 會一致）
    counts = accumarray(layer_of_orb, 1);
    if max(counts)-min(counts) > 0
        % 不強制報錯，僅提醒；仍可工作
        warning('各層無自旋軌域數不一致（min=%d, max=%d）。', min(counts), max(counts));
    end
    NorbL = mode(counts);

    % 3) 取 top/bot 的 Ns 層索引
    [~, order] = sort(z_levels, 'ascend');   % 由底到頂
    bot_layers = order(1:Ns);
    top_layers = order(end-Ns+1:end);

    mask_top_orb = ismember(layer_of_orb, top_layers);
    mask_bot_orb = ismember(layer_of_orb, bot_layers);

    % 4) 做投影器（先無自旋，再擴到含自旋）
    P_orb_top = spdiags(double(mask_top_orb), 0, Norb0, Norb0);
    P_orb_bot = spdiags(double(mask_bot_orb), 0, Norb0, Norb0);
    P_top     = kron(speye(2), P_orb_top);
    P_bot     = kron(speye(2), P_orb_bot);

    % info
    info.z = z; info.z_levels = z_levels; info.layer_of_orb = layer_of_orb;
    info.NL = NL; info.NorbL = NorbL; info.tol = tol;
end
