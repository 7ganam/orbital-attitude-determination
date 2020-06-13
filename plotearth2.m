function []=plotearth2(x,y,z,rad,im)
% imgRGB = imread('we5.png');
imgRGB = im;
[imgInd,map] = rgb2ind(imgRGB,256);
[imgIndRows,imgIndCols] = size(imgInd);
[X,Y,Z] = sphere(imgIndRows,imgIndCols);
surface((X*rad)+x,(Y*rad)+y,(Z*rad)+z,flipud(imgInd),...
    'FaceColor','texturemap',...
    'EdgeColor','none',...
    'CDataMapping','direct')
colormap(map)
% view(-35,45)
 axis vis3d
 text(2,2,2,'allp','HorizontalAlignment','left','FontSize',8);