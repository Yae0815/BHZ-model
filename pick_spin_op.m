function Sg = pick_spin_op(gamma, Norb)
    assert(mod(Norb,2)==0, 'Norb 必須為偶數');
    Norb_orb = Norb/2;
    switch lower(gamma)
        case 'x', s = [0 1; 1 0]/2;        % σ_x/2
        case 'y', s = [0 -1i; 1i 0]/2;     % σ_y/2
        case 'z', s = [1 0; 0 -1]/2;       % σ_z/2
        otherwise, error('gamma must be x/y/z');
    end
    Sg = kron(s, speye(Norb_orb));
end
