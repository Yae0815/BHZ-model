function write_spectral_weight_data(kx,nk,ky,EF,A,filename)
% kpath_car =linspace(-1,1,nk);
% A = sum(Ab_part(:,:,index),3);
% filename = [material_name '_' orbital_name '.dat'];

%% bulk or surface %%
fod = fopen(filename,'w');
fprintf(fod,' # of kx kpoints : %d \r\n',nk);
fprintf(fod,' kx_min : %11.6f ',min(kx));
fprintf(fod,' kx_max : %11.6f \r\n',max(kx));
fprintf(fod,' %11.6f ',kx);
fprintf(fod,'\r\n');
fprintf(fod,'\r\n');
fprintf(fod,' # of ky kpoints : %d \r\n',size(ky,2));
fprintf(fod,' ky_min : %11.6f ',min(ky-EF));
fprintf(fod,' ky_max : %11.6f \r\n',max(ky-EF));
fprintf(fod,' %11.6f ',ky-EF);
fprintf(fod,'\r\n');
fprintf(fod,'\r\n');
fprintf(fod,' spectral weight : %d * %d (row : kx, column : ky)\r\n',nk,size(ky,2));
for i_k = 1:nk
    fprintf(fod,' %11.6f ',A(:,i_k)/pi);
    fprintf(fod,'\r\n');
end
fclose(fod);