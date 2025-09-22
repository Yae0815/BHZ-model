function out = spin_hall_main(params)
% shc.spin_hall_main
% Compute spin Hall conductivity σ^{s_γ}_{αβ} using Kubo (AB−BA antisymmetrized).
% Supports T=0 step occupancy or finite-T Fermi function.
% Optional Kubo–Bastin weighting (f_n - f_m) for extra rigor.

    % --- unpack ---
    ftn   = params.ftn58;
    Ef    = params.Ef;
    Nk    = params.Nk;
    eta   = params.eta;
    e     = params.electronic_charge;
    hbar  = params.hbar;
    alpha = lower(params.alpha);
    beta  = lower(params.beta);
    gamma = lower(params.gamma);

    T   = 0;
    if isfield(params,'temperature_K'), T = params.temperature_K; end
    method = 'weighted';
    if isfield(params,'method'), method = lower(params.method); end
    kB = 8.617333262e-5; % eV/K

    % --- builders / sizes ---
    build = shc.make_builders(ftn);
    Norb  = build.Norb;
    Sg    = loc_pick_spin_op(gamma, Norb);   % 本地函式

    % --- Γ sampling (diagnostic) ---
    H_G = build.H(0,0,0);
    comm_SgH = Sg*H_G - H_G*Sg;
    comm_norm = norm(full(comm_SgH), 'fro');
    [Vg, Dg] = eig(full(H_G));
    Eg = real(diag(Dg)); [Eg,ix] = sort(Eg,'ascend'); Vg = Vg(:,ix);
    Sz_exp = real(diag(Vg'*(Sg*Vg)));
    sample_Gamma = [Eg, Sz_exp];

    % --- k-mesh (support shift) ---
    base = (0:Nk-1)/Nk - 0.5;
    shift = [0 0 0];
    if isfield(params,'shift'), shift = params.shift; end
    kx_list = base + shift(1)/Nk;
    ky_list = base + shift(2)/Nk;
    kz_list = base + shift(3)/Nk;
    w = 1.0 / (Nk^3);

    % --- capture handles ---
    H_of_k = build.H;
    dHx    = build.dHdkx;
    dHy    = build.dHdky;
    dHz    = build.dHdkz;

    % ===== parfor over kx slices =====
    sigma_partial = zeros(Nk,1);

    parfor ixk = 1:Nk
        kx = kx_list(ixk);
        local_sum = 0.0;

        for iyk = 1:Nk
            ky = ky_list(iyk);
            for izk = 1:Nk
                kz = kz_list(izk);

                % H and velocities
                Hk  = H_of_k(kx,ky,kz);
                vA  = loc_pick_velocity(dHx,dHy,dHz, alpha, kx,ky,kz, hbar);
                vB  = loc_pick_velocity(dHx,dHy,dHz, beta,  kx,ky,kz, hbar);

                % numerical hermitization
                vA  = (vA + vA')/2;   vB  = (vB + vB')/2;
                JAs = 0.5*(Sg*vA + vA*Sg);
                JBs = 0.5*(Sg*vB + vB*Sg);
                JAs = (JAs + JAs')/2; JBs = (JBs + JBs')/2;

                % diag once
                [U,D] = eig(full(Hk));
                E = real(diag(D)); [E,idx] = sort(E,'ascend'); U = U(:,idx);

                % eigenbasis
                J_A = U'*(JAs*U);  J_B = U'*(JBs*U);
                v_A = U'*(vA *U);  v_B = U'*(vB *U);

                % Kubo, antisymmetrized
                Om_AB = 0.0;  % parfor needs explicit init
                Om_BA = 0.0;

                switch method
                    case 'weighted' % Fermi-sea with Fermi weights w_n
                        if T<=0
                            w_n = double(E < Ef);
                        else
                            w_n = 1 ./ (1 + exp((E - Ef)/(kB*T)));
                        end
                        Om_AB = loc_pair_weighted(J_A, v_B, E, eta, w_n);
                        Om_BA = loc_pair_weighted(J_B, v_A, E, eta, w_n);

                    case 'bastin'   % rigorous (f_n - f_m) weighting
                        if T<=0
                            Om_AB = loc_pair_bastin(J_A, v_B, E, eta, Ef, 0);
                            Om_BA = loc_pair_bastin(J_B, v_A, E, eta, Ef, 0);
                        else
                            Om_AB = loc_pair_bastin(J_A, v_B, E, eta, Ef, T);
                            Om_BA = loc_pair_bastin(J_B, v_A, E, eta, Ef, T);
                        end

                    otherwise
                        error('params.method must be ''weighted'' or ''bastin''.');
                end

                Om = 0.5*(Om_AB - Om_BA);
                local_sum = local_sum + Om * w;
            end
        end

        sigma_partial(ixk) = local_sum;
    end

    sigma_val = (e/hbar) * sum(sigma_partial);

    % --- outputs ---
    out.sigma_sab_gamma = sigma_val;
    out.alpha = alpha; out.beta = beta; out.gamma = gamma;
    out.commutator_norm_Gamma = comm_norm;
    out.sample_bands_Gamma = sample_Gamma;
    out.note = 'Reduced-k lattice units; prefactor e/ħ applied.';
end

% ==== Local helpers (no path/look-up headaches) ===========================

function Sg = loc_pick_spin_op(gamma, Norb)
    assert(mod(Norb,2)==0, 'Norb must be even (spin×orbital).');
    Norb_orb = Norb/2;
    switch lower(gamma)
        case 'x', s = [0 1; 1 0]/2;
        case 'y', s = [0 -1i; 1i 0]/2;
        case 'z', s = [1 0; 0 -1]/2;
        otherwise, error('gamma must be x/y/z');
    end
    Sg = kron(s, speye(Norb_orb));
end

function v = loc_pick_velocity(dHx,dHy,dHz, alpha, kx, ky, kz, hbar)
    twoPi = 2*pi;
    switch lower(alpha)
        case 'x', dH = dHx(kx,ky,kz);
        case 'y', dH = dHy(kx,ky,kz);
        case 'z', dH = dHz(kx,ky,kz);
        otherwise, error('alpha must be x/y/z');
    end
    v = (1/hbar) * (1/twoPi) * dH;
end

function Om = loc_pair_weighted(Jmat, vmat, E, eta, w_n)
    Norb = numel(E);
    Om = 0.0;
    for n = 1:Norb
        if w_n(n)==0, continue; end
        dE2 = (E - E(n)).^2 + eta^2;
        dE2 = max(dE2, eta^2);
        num = Jmat(n,:).*vmat(:,n).';
        num(n) = 0;
        dE2(n) = Inf;
        Om = Om + w_n(n) * 2*imag( sum( num ./ dE2.' ) );
    end
end

function Om = loc_pair_bastin(Jmat, vmat, E, eta, mu, T)
    kB = 8.617333262e-5; % eV/K
    if T<=0
        f = @(x) double(x < mu);
    else
        f = @(x) 1./(1+exp((x-mu)/(kB*T)));
    end
    Norb = numel(E);
    Om = 0.0;
    for n = 1:Norb
        fn = f(E(n));
        for m = 1:Norb
            if m==n, continue; end
            fm = f(E(m));
            denom = (E(m)-E(n))^2 + eta^2;
            Om = Om + (fn - fm) * 2*imag( Jmat(n,m)*vmat(m,n) ) / denom;
        end
    end
end
