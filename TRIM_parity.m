%% TRIM Parity Calculator for ftn58sparse Hamiltonian
%   ---------------------------------------------------
%   * Computes inversion parity at the 8 TRIM and Fu–Kane ℤ₂ indices.
%   * δ_Λ uses **one state per Kramers pair**, selected by taking **odd
%     band indices (1,3,5,…) among occupied levels** – a basis‑agnostic
%     rule that remains valid as long as the spin‑degenerate partners are
%     stored consecutively (↑ then ↓).
%   ---------------------------------------------------
%   Usage
%   -----
%   1. Place with **ftn58sparse.mat**.
%   2. Run:  >> trim_parity
%   3. Results appear on screen and in **trim_parity.mat**.
%   ---------------------------------------------------
%   Last update: Aug‑2025 (δ_Λ selects odd occupied band indices)
%   ---------------------------------------------------

clear;  clc;  tic;

%% (0)  Parameters  ----------------------------------------------------
Ef = 0;                         % Fermi level (eV)

%% (1)  Load tight‑binding data  --------------------------------------
load ftn58sparse.mat            % norb, ij, dd, tt, …

norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
dd   = ftn58sparse.dd;
tt   = ftn58sparse.tt;

%% (2)  Inversion operator P  -----------------------------------------
Porb = diag([ 1  -1   1  -1]);          % orbital‑level parity
Pbig = kron(Porb, eye(norb/4));         % full operator

%% (3)  8 TRIM  --------------------------------------------------------
trim_frac = [ ...
    0   0   0 ;   0.5 0   0 ;   0   0.5 0 ;   0   0   0.5 ;
    0.5 0.5 0 ;   0.5 0   0.5 ; 0   0.5 0.5 ; 0.5 0.5 0.5 ];
trim_lab  = { 'Γ' 'X' 'Y' 'Z' 'S' 'U' 'T' 'R' };
kpts      = trim_frac * 2*pi;

%% (4)  Storage  -------------------------------------------------------
Eigs_all    = zeros(norb,8);
Parity_all  = zeros(norb,8,'int8');
Delta_TRIM  = zeros(1,8,'int8');

%% (5)  Loop over TRIM  ------------------------------------------------
for ik = 1:8
    k = kpts(ik,:).';
    H = sparse(ii, jj, exp(1i*dd*k).*tt, norb, norb);
    H = (full(H)+full(H)')/2;           % Hermitian safety
    [V,D] = eig(H,'vector');

    Eigs_all(:,ik)   = real(D);
    Parity_all(:,ik) = int8(sign(real(diag(V' * Pbig * V))));

    occ_mask = (D <= Ef + 1e-6);        % Fermi‑filled bands
    band_idx = (1:norb).';
    odd_idx  = mod(band_idx,2)==1;      % 1,3,5,… — one per Kramers pair
    sel_mask = occ_mask & odd_idx;

    Delta_TRIM(ik) = prod(Parity_all(sel_mask,ik));
end

toc;

%% (6)  Print δ_Λ  -----------------------------------------------------
fprintf('\n========  TRIM Parity Summary  ========\n');
for ik = 1:8
    fprintf('  %s (δ_%s = %2d)\n', trim_lab{ik}, trim_lab{ik}, Delta_TRIM(ik));
end

%% (7)  Fu–Kane Z₂  ----------------------------------------------------
nu0 = mod(sum(Delta_TRIM == -1),2);
nu1 = mod(sum(Delta_TRIM(trim_frac(:,1)==0.5) == -1),2);
nu2 = mod(sum(Delta_TRIM(trim_frac(:,2)==0.5) == -1),2);
nu3 = mod(sum(Delta_TRIM(trim_frac(:,3)==0.5) == -1),2);

fprintf('\nFu–Kane ℤ₂  →  (ν₀; ν₁ ν₂ ν₃) = (%d; %d %d %d)\n', nu0, nu1, nu2, nu3);

%% (8)  Save  ----------------------------------------------------------
results.trim_frac  = trim_frac;
results.kpts       = kpts;
results.Ef         = Ef;
results.Eigs_all   = Eigs_all;
results.Parity_all = Parity_all;
results.Delta_TRIM = Delta_TRIM;
results.Z2         = [nu0 nu1 nu2 nu3];

save('trim_parity.mat','results');

fprintf('\nData saved: trim_parity.mat\n');