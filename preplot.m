function[earth]=preplot(axes3,a,b)
%----------------------------
a;
b;
axis equal
rad = 6378.136;
figure(3)       
axes(axes3)

% plot axis-----------------------------------------
hold on
h1 =  plot3([0 rad+5000],[0 0],[0 0] ,'r-');
h2 =  plot3([0 0],[0 rad+5000],[0 0] ,'g-');
h3 =  plot3([0 0],[0 0],[0 rad+5000],'color',[0 0 .8]);
set([h1 h2 h3],'linewidth',4);
view(45,20);
text(rad+6000,0,0,'x','HorizontalAlignment','left','FontSize',10);
text(0,rad+6000,0,'y','HorizontalAlignment','left','FontSize',10);
text(0,0,rad+6000,'z','HorizontalAlignment','left','FontSize',10);
axis([-a a -b b])
% camva(5);
  
%--------------------------------------------------

% plot equtorial plan -------------------------------
yu=100000;
pointa1 = [-2*yu,-2*yu,0];
pointa2 = [-2*yu,2*yu,0];
pointb1 = [2*yu,2*yu,0];
pointb2 = [2*yu,-2*yu,0];
pointbs=[pointa1' pointa2' pointb1' pointb2']; % using the data given in the question
fill3(pointbs(1,:),pointbs(2,:),pointbs(3,:),'b')
alpha(0.5)

ang=0:0.01:2*pi; 
xp=rad*cos(ang);
yp=rad*sin(ang);
plot(xp,yp,'b');
% ---------------------------------------------------

% plot earth --------------------------------------6
earth=plotearth(0,0,0,rad);
direction = [0 0 1];
rotate(earth,direction,0)
axis equal
% -------------------------------------------------
