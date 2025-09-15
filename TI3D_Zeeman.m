%%% Program to Generate ftn58 Fromat from Lattice Tight Binding Model For Gamma matrix TI model %%%
%%% Reference: Rundong Li et al., Nat. Phys. 6, 284-288 (2010)                                  %%%
%%% Author: Yueh-Ting Yao (2024/01/05)                                                          %%%

clear all

%%-- Parameters --%%
A1 = 0.4000;  A2 = 0.4000;
B1 = 0.3;  B2x = 0.3;  B2y = 0.3;
C = 0.0000;  D1 = 0.0000;  D2 = 0.0000;  M = 1;
  DO = 0.45;  Delta = 0.00;  Rx = 0.00; Ry = 0.00;

% ----------------------------------------------------------------------------------
%%-- Generating Hopping Terms --%%
%% onsite
ftn58_onsite_11 = [[1 1 (C+2*D1+2*D2) 0 0 0];[1 1 -D1 0 0 1];[1 1 -D2 1 0 0];[1 1 -D2 0 1 0];...
                                             [1 1 -D1 0 0 -1];[1 1 -D2 -1 0 0];[1 1 -D2 0 -1 0]];
ftn58_onsite_22 = ftn58_onsite_11;
ftn58_onsite_22(:,1:2) = ftn58_onsite_11(:,1:2) + 1;
ftn58_onsite_33 = ftn58_onsite_11;
ftn58_onsite_33(:,1:2) = ftn58_onsite_11(:,1:2) + 2;
ftn58_onsite_44 = ftn58_onsite_11;
ftn58_onsite_44(:,1:2) = ftn58_onsite_11(:,1:2) + 3;
ftn58_onsite    = [ftn58_onsite_11;ftn58_onsite_22;ftn58_onsite_33;ftn58_onsite_44];

%% Gamma 1
ftn58_G1 = [[1 4 A2*(1/2i) 1 0 0];[2 3 A2*(1/2i) 1 0 0];...
            [1 4 (-1)*A2*(1/2i) -1 0 0];[2 3 (-1)*A2*(1/2i) -1 0 0]];

%% Gamma 2
ftn58_G2 = [[1 4 (-1i)*A2*(1/2i) 0 1 0];[2 3 (-1i)*A2*(1/2i) 0 1 0];...
            [1 4 (-1)*(-1i)*A2*(1/2i) 0 -1 0];[2 3 (-1)*(-1i)*A2*(1/2i) 0 -1 0]];

%% Gamma 3
ftn58_G3 = [[1 2 (-1i)*A1*(1/2i) 0 0 1];[3 4 (-1i)*A1*(1/2i) 0 0 1];...
            [1 2 (-1)*(-1i)*A1*(1/2i) 0 0 -1];[3 4 (-1)*(-1i)*A1*(1/2i) 0 0 -1]];

%% Gamma 4
ftn58_G4_11 = [[1 1 (M-2*B1-2*B2x-2*B2y) 0 0 0];[1 1 B1 0 0 1];[1 1 B2x 1 0 0];[1 1 B2y 0 1 0];...
                                         [1 1 B1 0 0 -1];[1 1 B2x -1 0 0];[1 1 B2y 0 -1 0]];
ftn58_G4_22 = ftn58_G4_11;
ftn58_G4_22(:,1:2) = ftn58_G4_11(:,1:2) + 1;
ftn58_G4_22(:,3)   = (-1)*ftn58_G4_11(:,3);
ftn58_G4_33 = ftn58_G4_11;
ftn58_G4_33(:,1:2) = ftn58_G4_11(:,1:2) + 2;
ftn58_G4_44 = ftn58_G4_11;
ftn58_G4_44(:,1:2) = ftn58_G4_11(:,1:2) + 3;
ftn58_G4_44(:,3)   = (-1)*ftn58_G4_11(:,3);
ftn58_G4    = [ftn58_G4_11;ftn58_G4_22;ftn58_G4_33;ftn58_G4_44];

%% Zeeman
ftn58_zm = [[1 1 Delta 0 0 0];[2 2 Delta 0 0 0];[3 3 -Delta 0 0 0];[4 4 -Delta 0 0 0]];

%% Rashba
ftn58_Rashba = [[1 3 (1/2)*Ry 0 1 0];[1 3 (-1/2)*Ry 0 -1 0];[1 3 1i*(1/2)*Rx 1 0 0];[1 3 1i*(-1/2)*Rx -1 0 0];[2 4 (1/2)*Ry 0 1 0];[2 4 (-1/2)*Ry 0 -1 0];[2 4 1i*(1/2)*Rx 1 0 0];[2 4 1i*(-1/2)*Rx -1 0 0];];


%% inter-orbit coupling
ftn58_DO = [[1 2 DO 0 0 0];[3 4 DO 0 0 0]];

% ----------------------------------------------------------------------------------

ftn58 = [ftn58_onsite;ftn58_G1;ftn58_G2;ftn58_G3;ftn58_G4;ftn58_zm;ftn58_Rashba;ftn58_DO];
nbond = length(ftn58);
ftn58 = [(1:nbond)' ftn58];
ftn58 = [[4 nbond 0 0 0 0 0 ];ftn58];

             
save ftn58.mat ftn58
eval('run TBHmftn.m');
eval('run band_ftn.m');

