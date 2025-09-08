function [SxP, SyP, SzP, W] = projected_spin(vec, Sx_full, Sy_full, Sz_full, P)
% 對每個本徵態做「投影後歸一化」的自旋期望值與權重
%   vec:  (norb x norb) 本徵向量（欄向量逐一為本徵態）
%   S?_full: 含自旋的自旋矩陣（例如 kron(Sx, eye(norb/2))）
%   P:    (norb x norb) 投影器（上或下表面）
% 回傳：
%   SxP,SyP,SzP: (1 x norb) 每條能帶投影自旋
%   W:           (1 x norb) 投影權重（該態落在投影子空間的機率）
    % 權重
    WT = vec' * (P*vec);   WT = 0.5*(WT+WT');  w = real(diag(WT)).';
    % 各向自旋（投影後）
    num = vec' * (P*Sx_full*P) * vec; num = 0.5*(num+num'); SxP = real(diag(num)).';
    num = vec' * (P*Sy_full*P) * vec; num = 0.5*(num+num'); SyP = real(diag(num)).';
    num = vec' * (P*Sz_full*P) * vec; num = 0.5*(num+num'); SzP = real(diag(num)).';
    % 歸一化（避免 0 除）
    denom = max(w, 1e-12);
    SxP = SxP ./ denom;  SyP = SyP ./ denom;  SzP = SzP ./ denom;
    W   = w;
end
