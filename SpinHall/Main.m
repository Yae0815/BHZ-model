%% ===== Main.m : compute full spin Hall tensor σ^{s_γ}_{αβ} =====
clear; clc;

% --- Load TB model ---
load ../ftn58sparse.mat

% --- Base params ---
params.ftn58 = ftn58sparse;
params.Ef    = 0.0;
params.Nk    = 41;       % odd
params.eta   = 1e-4;
params.electronic_charge = 1.0;
params.hbar  = 1.0;

axes = 'xyz';
tensor = zeros(3,3,3);      % gamma x alpha x beta
commGamma = zeros(3,1);

fprintf('Running unique σ^{s_γ}_{αβ} (α<β) ... Nk=%d, eta=%.1e\n', params.Nk, params.eta);
t0 = tic;

for gi = 1:3
    params.gamma = axes(gi);
    got_comm = false;

    for ai = 1:3
        params.alpha = axes(ai);

        % 只算 α<β
        for bi = ai+1:3
            params.beta = axes(bi);

            res = spin_hall_main(params);  % 內部已 parfor + 反對稱化
            val = res.sigma_sab_gamma;

            % 寫入上三角，並用反對稱性補齊下三角
            tensor(gi, ai, bi) =  val;
            tensor(gi, bi, ai) = -val;

            fprintf('sigma^{s_%s}_{%s%s} = % .10e\n', params.gamma, params.alpha, params.beta, val);

            if ~got_comm
                commGamma(gi) = res.commutator_norm_Gamma;
                got_comm = true;
            end
        end
    end
end

elapsed = toc(t0);
fprintf('\nDone. Elapsed time: %.2f s\n', elapsed);

% --- Print 3 matrices (one per γ), and antisymmetry checks ---
labels = {'x','y','z'};
for gi = 1:3
    M = squeeze(tensor(gi,:,:));
    M(1:4:9) = 0;
    fprintf('\n=== σ^{s_%s}_{αβ} (3×3) ===\n', axes(gi));
    print_matrix_with_labels(M, labels);
    as_err = norm(M + M.', 'fro');     % 應接近 0（數值誤差）
    fprintf('antisymmetry check  ||M + M^T||_F = %.3e\n', as_err);
    fprintf('||[S_%s, H(Γ)]||_F = %.3e\n', axes(gi), commGamma(gi));
end

% --- Unique components table (α<β) ---
unique_list = {'xy','xz','yz'};
U = zeros(3,3);
U(:,1) = [tensor(1,1,2); tensor(2,1,2); tensor(3,1,2)];
U(:,2) = [tensor(1,1,3); tensor(2,1,3); tensor(3,1,3)];
U(:,3) = [tensor(1,2,3); tensor(2,2,3); tensor(3,2,3)];
fprintf('\n=== Unique components (α<β): rows γ∈{x,y,z}, cols {xy, xz, yz} ===\n');
print_matrix_with_labels(U, unique_list);

% --- Save results ---
meta = struct('Ef',params.Ef,'Nk',params.Nk,'eta',params.eta, ...
              'units','lattice (e=ħ=1); reduced k; BZ vol=1', ...
              'elapsed_sec',elapsed);
save('spin_hall_tensor.mat','tensor','meta');
fprintf('\nSaved: spin_hall_tensor.mat\n');

% helper
function print_matrix_with_labels(M, labels)
    if ischar(labels), labels = cellstr(labels.'); end
    fprintf('        %6s         %6s         %6s\n', labels{1}, labels{2}, labels{3});
    for i = 1:3
        fprintf('  %s | % .10e  % .10e  % .10e\n', labels{i}, M(i,1), M(i,2), M(i,3));
    end
end
