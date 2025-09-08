function v = pick_velocity(build, alpha, kx,ky,kz, hbar)
    twoPi = 2*pi;
    switch lower(alpha)
        case 'x', dH = build.dHdkx(kx,ky,kz);
        case 'y', dH = build.dHdky(kx,ky,kz);
        case 'z', dH = build.dHdkz(kx,ky,kz);
        otherwise, error('alpha must be x/y/z');
    end
    v = (1/hbar) * (1/twoPi) * dH;
end