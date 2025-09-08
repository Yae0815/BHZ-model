function cm_data=blue_red(m)
colormap_00 = [
    255,   0,   0;
    229 53 11;
   255,   255,   255;
   0 176 240;
   0,   0, 255 ];
colormap_0 = colormap_00/255.0;
blue_red = interp1(1:size(colormap_0,1), ...
colormap_0, ...
linspace(1,size(colormap_0,1),64), 'pchip');
cm_data = blue_red;
end