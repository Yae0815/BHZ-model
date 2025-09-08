Sx = [0 1; 1 0];
Sy = [0 -1i; 1i 0];
Sz = [1 0; 0 -1];
S0 = eye(2);
Gamma_3 = kron(S0, Sy);
Gamma_5 = kron(Sz, Sx);

Gamma_35 = (1/2i)*(Gamma_3*Gamma_5-Gamma_5*Gamma_3);

Gamma_4 = kron(S0,Sz);

%disp(Gamma_35);
disp(Gamma_35*Gamma_4-Gamma_4*Gamma_35);