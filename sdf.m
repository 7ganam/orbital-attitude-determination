for t=0:.01:2*pi
 X = a*cos(t);
 Y = b*sin(t);
 w = atan2(0,2*c);
 x =  X*cos(w) - Y*sin(w);
 y =  X*sin(w) + Y*cos(w);
plot(x,y,'-.r*');
w=r/10;
imagesc([c-r,c+r],[-r,r],im)
pause(.01)
 drawnow
end