% demo_mu_diagram.m  (precompute + μ-scan)
clear; clc;

% ===== add path =====
this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(this_dir, '..'))); % add shc_kubo/+shc

% ===== load model =====
ftn = load('../ftn58sparse.mat');
if isfield(ftn,'ftn58sparse'), ftn = ftn.ftn58sparse; end

% ===== params =====
params.ftn58              = ftn;
params.Nk                 = 41;
params.eta                = 5e-3;           % eV (fixed in cache; change => re-precompute)
params.electronic_charge  = 1;              % lattice units
params.hbar               = 1;              % lattice units
params.alpha              = 'x';
params.beta               = 'z';
params.gamma              = 'x';
params.shift              = [0 0 0];        % can test [0.5 0.5 0.5]/Nk
% (Ef不重要，評分時會給 mu)

% ===== 1) precompute once =====
cache = shc.precompute_kgrid(params);

% ===== 2) evaluate for many μ (and T if wanted) =====
mu_grid = linspace(-0.5, 0.5, 101);
T       = 0;                  % change to e.g. 300 (K) to test finite-T
method  = 'weighted';         % or 'bastin'

sig = zeros(numel(mu_grid),1);
for i = 1:numel(mu_grid)
    out = shc.eval_sigma(cache, mu_grid(i), T, method);
    sig(i) = out.sigma;
end

% ===== 3) plot =====
figure('Color','w'); plot(mu_grid, sig, 'LineWidth',1.6);

% labels (純 lattice units，去掉 e/ħ)
lbl_interpreter = 'latex';                      % 'latex' 或 'tex'
unit_note       = '\mathrm{(lattice\ units)}';

a = lower(params.alpha); b = lower(params.beta); g = lower(params.gamma);
switch lower(lbl_interpreter)
    case 'latex'
        yl_core = sprintf('\\sigma^{s_{%s}}_{%s%s}', g, a, b);
        yl = ['$' yl_core '\; ' unit_note '$'];
        xl = '$\mu\ \mathrm{(eV)}$';
    otherwise
        yl = sprintf('\\sigma^{s_{%s}}_{%s%s} %s', g, a, b, unit_note);
        xl = '\mu (eV)';
end
xlabel(xl, 'Interpreter', lbl_interpreter, 'FontSize', 14);
ylabel(yl, 'Interpreter', lbl_interpreter, 'FontSize', 14);
title('Spin Hall vs. Chemical Potential', 'Interpreter', 'none');
grid on; box on;

fprintf('Precompute done: Nk=%d, eta=%.3g eV. Eval method=%s, T=%.1f K\n', ...
    params.Nk, params.eta, method, T);
