% @author Angus Huang rabbit4a9@gmail.com
function ink_blue = colormap_ink_blue


colormap_00 = [
    255, 255, 255;
     77,  77, 244;
      0,   0, 202;
      0,   0, 111;
      0,   0,   0];


colormap_0 = colormap_00/255.0;
ink_blue = interp1(1:size(colormap_0,1), ...
    colormap_0, ...
    linspace(1,size(colormap_0,1),1024), 'pchip');

