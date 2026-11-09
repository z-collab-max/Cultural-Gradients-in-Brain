mycolorpoint = [
    [103,0,31];
    [178,24,43];
    [214,96,77];
    [244,165,130];
    [253,219,199];
    [247,247,247];
    [209,229,240];
    [146,197,222];
    [67,147,195];
    [33,102,172];
    [5,48,97]];
mycolorposition=[1    12    23    34    45    56    67    78    89   100   111];
mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),1:111,'linear','extrap');
mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),1:111,'linear','extrap');
mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),1:111,'linear','extrap');
mycolor=[mycolormap_r',mycolormap_g',mycolormap_b']/255;
mycolor=round(mycolor*10^4)/10^4;
mycolor = flip(mycolor);
