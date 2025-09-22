function build = make_builders(ftn)
% shc.make_builders
% Return H(k) and ∂H/∂kα(k) handles for reduced k (phase e^{i 2π k·R}).

    ij = ftn.ij;   % [nbond x 2]
    tt = ftn.tt;   % [nbond x 1] complex
    dd = ftn.dd;   % [nbond x 3] integers
    Norb = ftn.norb;

    i_list = ij(:,1);
    j_list = ij(:,2);
    dx = dd(:,1); dy = dd(:,2); dz = dd(:,3);

    twoPi = 2*pi;

    function Hk = H_of_k(kx,ky,kz)
        phase = exp(1i*twoPi*(kx*dx + ky*dy + kz*dz));
        val = tt .* phase;
        Hk = accumSparseHermitian(i_list, j_list, val, Norb);
    end

    function dHx = dHdkx(kx,ky,kz)
        phase = exp(1i*twoPi*(kx*dx + ky*dy + kz*dz));
        val = (1i*twoPi*dx) .* tt .* phase;
        dHx = accumSparseHermitian(i_list, j_list, val, Norb);
    end

    function dHy = dHdky(kx,ky,kz)
        phase = exp(1i*twoPi*(kx*dx + ky*dy + kz*dz));
        val = (1i*twoPi*dy) .* tt .* phase;
        dHy = accumSparseHermitian(i_list, j_list, val, Norb);
    end

    function dHz = dHdkz(kx,ky,kz)
        phase = exp(1i*twoPi*(kx*dx + ky*dy + kz*dz));
        val = (1i*twoPi*dz) .* tt .* phase;
        dHz = accumSparseHermitian(i_list, j_list, val, Norb);
    end

    build.H     = @H_of_k;
    build.dHdkx = @dHdkx;
    build.dHdky = @dHdky;
    build.dHdkz = @dHdkz;
    build.Norb  = Norb;
end

function H = accumSparseHermitian(i_list, j_list, val, Norb)
    H = sparse(i_list, j_list, val, Norb, Norb);
    H = (H + H')/2; % numeric symmetrization
end
