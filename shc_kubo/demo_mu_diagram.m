% demo_mu_diagram.m
% Quick start for σ^{s_γ}_{αβ}(μ) diagram

clear; clc;

% ===== 0) add path (only the package root) =====
% Change this to your actual location:
this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(this_dir, '..'))); % add shc_kubo/+shc to path

% ===== 1) load model =====
% Expecting ../ftn58sparse.mat with fields: ij, tt, dd, norb
% (You can change the relative path below)
ftn = load('../ftn58sparse.mat');
if isfield(ftn,'ftn58sparse')
    ftn = ftn.ftn58sparse;
end

% ===== 2) parameters =====
params.ftn58              = ftn;
params.Ef                 = 0.0;            % will be overwritten per μ
params.Nk                 = 41;             % increase for convergence
params.eta                = 5e-3;           % eV
%params.electronic_charge  = 1.602176634e-19; % C
params.electronic_charge  = 1;
%params.hbar               = 6.582119569e-16; % eV*s
params.hbar               = 1;
params.alpha              = 'x';
params.beta               = 'z';
params.gamma              = 'x';
params.shift              = [0 0 0];        % try [0.5 0.5 0.5] / Nk for smoothing
params.temperature_K      = 0;              % set >0 for finite T
params.method             = 'weighted';     % 'weighted' or 'bastin'

% ===== 3) μ grid =====
mu_grid = linspace(-0.5, 0.5, 3);         % eV; adapt to your band range

% ===== 4) run scan (no nested parfor inside) =====
scan = shc.scan_mu(params, mu_grid);

% ===== 5) plot =====
figure('Color','w');
plot(scan.mu, scan.sigma, 'LineWidth',1.6);

% 可調整的標籤參數
lbl_interpreter = 'latex';                      % 'latex' 或 'tex'
unit_note       = '\mathrm{(lattice\ units)}';  % 單位註記（LaTeX 下底線已轉義）

% 依 params.alpha/beta/gamma 自動組字串
a = lower(params.alpha); b = lower(params.beta); g = lower(params.gamma);

switch lower(lbl_interpreter)
    case 'latex'
        % 整段放進數學模式，避免 LaTeX 指令在文字模式報錯
        yl_core = sprintf('\\sigma^{s_{%s}}_{%s%s}', g, a, b);
        yl = ['$' yl_core '\; (e/\hbar \cdot ' unit_note ')$'];
        xl = '$\mu\ \mathrm{(eV)}$';   % 參考你的範例：ylabel('$k_y\,(2\pi/a)$',...)
    case 'tex'
        % TeX 模式不需 $...$，也不支援 \mathrm 與 \cdot（可留空白代替）
        yl = sprintf('\\sigma^{s_{%s}}_{%s%s} (e/\\hbar  %s)', g, a, b, unit_note);
        xl = '\mu (eV)';
    otherwise
        error('lbl_interpreter must be ''latex'' or ''tex''.');
end

xlabel(xl, 'Interpreter', lbl_interpreter, 'FontSize', 14);
ylabel(yl, 'Interpreter', lbl_interpreter, 'FontSize', 14);

title('Spin Hall vs. Chemical Potential', 'Interpreter', 'none');
grid on; box on;

% ===== 6) quick prints =====
fprintf('Done. Nk=%d, eta=%.3g eV, T=%.1f K, method=%s\n', ...
    params.Nk, params.eta, params.temperature_K, params.method);
