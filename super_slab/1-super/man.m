clear
    filename = ['S2F' ,num2str(ii), '/ftn58soc_full.mat'];
    filenameT = ['S2F' ,num2str(ii), '/G_dist.mat'];
    eval(['load ' filename]);
    ftn58=ftn58soc_full;
%% Kpoints %%%
%    nk = 100;
%    p1 = [linspace(0.0,1/3,nk+1)' linspace(0,1/3,nk+1)'];
%    p2 = [linspace(1/3,0.5,nk+1)' linspace(1/3,0.0,nk+1)'];
%    p3 = [linspace(0.5,0.0,nk)' linspace(0.0,0.0,nk)'];

%    list    = [1 101 201 300];
%    kpoints = [p1(1:nk,:);p2(1:nk,:);p3(1:nk,:)];
%    kpoints = [kpoints zeros(length(kpoints(:,1)),1)];
%% Gap check %%
    nk = 200;
    p1 = [linspace(-0.1,0,nk+1)' linspace(0,0,nk+1)'];
    p2 = [linspace(0,0.1,nk)' linspace(0,0,nk)'];

    list    = [1 201 400];
    kpoints = [p1(1:nk,:);p2(1:nk,:)];
    kpoints = [kpoints zeros(length(kpoints(:,1)),1)];
    
    [A,vecnorm]=bandplot(ftn58,kpoints);
    eval(['save ' filenameT ' A vecnorm']);
