function cm_data=red(m)
colormap_00 = [
    255,   255,   255;
    229 53 11;
    255,   0,   0];
colormap_0 = colormap_00/255.0;
blue_red = interp1(1:size(colormap_0,1), ...
colormap_0, ...
linspace(1,size(colormap_0,1),64), 'pchip');
cm_data = blue_red;
end