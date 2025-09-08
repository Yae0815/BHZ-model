for No=8
    for ii = 20
        clearvars -except ii No
        if No ==1
            filenameNorb = 'Cr-dz2';
        elseif No ==2
            filenameNorb = 'Cr-dxz,Cr-dyz';
        elseif No ==3
            filenameNorb = 'Cr-d(x2-y2),Cr-dxy';
        elseif No ==4
            filenameNorb = 'Te-pz';
        elseif No ==5
            filenameNorb = 'Te-px,Te-py';
        elseif No ==6
            filenameNorb = 'up';
        elseif No ==7
            filenameNorb = 'dn';
        elseif No ==8
            filenameNorb = 'norm';
        end
        %         filename = ['S2F' ,num2str(ii), '/ftn58soc_full.mat'];
        filenameT = ['S2F' ,num2str(ii), '/G_dist.mat'];
        filenameS = ['S2F' ,num2str(ii), '/G_' filenameNorb '.png'];
        filenameJ = ['S2F' ,num2str(ii), '/G_' filenameNorb '.fig'];
        eval(['load ' filenameT]);
        Ef = 5.15497639;
        TT=ii*24;
        %%
        if No == 1
            norb = [2,13];
        elseif No == 2
            norb = [3:4,14:15];
        elseif No == 3
            norb = [5:6,16:17];
        elseif No == 4
            norb = [7,10,19,22];
        elseif No == 5
            norb = [8:9,11:12,20:21,23:24];
        elseif No ==6
            norb=1:TT/2;
        elseif No ==7
            norb=(TT/2+1):TT;
        elseif No ==8
            norb = 1:TT;
        end
        %% Gap check %%
        nk = 200;
        p1 = [linspace(-0.1,0,nk+1)' linspace(0,0,nk+1)'];
        p2 = [linspace(0,0.1,nk)' linspace(0,0,nk)'];
        
        list    = [1 201 400];
        kpoints = [p1(1:nk,:);p2(1:nk,:)];
        kpoints = [kpoints zeros(length(kpoints(:,1)),1)];
        %% Kpoints-2
%         nk = 100;
%         p1 = [linspace(0.0,1/3,nk+1)' linspace(0,1/3,nk+1)'];
%         p2 = [linspace(1/3,0.5,nk+1)' linspace(1/3,0.0,nk+1)'];
%         p3 = [linspace(0.5,0.0,nk)' linspace(0.0,0.0,nk)'];
%         list    = [1 101 201 300];
%         kpoints = [p1(1:nk,:);p2(1:nk,:);p3(1:nk,:)];
        
        figure
        charA=vecnorm(norb,:,:);
        charA=sum(charA,1);
        hold on;
        nknb=size(A(:,1));
        if No ==8
            Ak = reshape(A(:,2)-Ef,length(norb),length(kpoints(:,1)));
            Ak = Ak';
            kd = 1:length(kpoints(:,1));
            plot(kd,Ak,'k-','LineWidth',1.5);
        else
            scatter(A(:,1),A(:,2)-Ef,10,reshape(charA,nknb(1),1),'filled');
            colorbar('vertical')
            % title(['orb:s']);
            caxis([0 1]);
            ym=ylim;
            xm=xlim;
            zm=caxis;
            plot(xm,[0 0],'-.');
            colormap(copper)
            box on
            axis square
        end
        %%
%         axis([0 300 -2 2]);
%         xticks([1 101 201 300]);
%         xm=xlim;
%         plot(xm,[0 0],'-.');
%         xticklabels({'\Gamma','K','M','\Gamma'});
%         ylabel('E-E_F(eV)')
%         set(gca,'FontSize',20)
        %%
        axis([0 400 1.4 1.6]);
        xticks([1 201 400]);
        xm=xlim;
        plot(xm,[0 0],'-.');
        xticklabels({'-M','\Gamma','+M'});
        ylabel('E-E_F(eV)')
        set(gca,'FontSize',20)
        
        saveas(gcf,filenameS);
        saveas(gcf,filenameJ);
        %         close all
    end
end