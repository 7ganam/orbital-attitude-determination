clear
 rotate3d off
 subplot(2,1,1)       
 axes1=gca;
subplot(2,1,2)       
 axes2=gca;
cla(axes1,'reset');
cla(axes2,'reset');

day=24*60*60;


% x=-6634;
x=7000;
y=7000;
z=2;
vx=-5;
vy=5;
vz=.5;

% assignement main code
% input V & R 
mu=3.9868*10^5;
V=[vx vy vz];
R= [x y z];
% example v=[  8 9 2];
% R=[ 7 2 9];

% find the magnitude of R and V
v=sqrt(V(1)^2+V(2)^2+V(3)^2);
r=sqrt(R(1)^2+R(2)^2+R(3)^2);

% clculating the specific energy ertyu
E=(v^2/2)-(mu/r);

% calculating specific angular momentu h
H= cross(R,V);
h=sqrt(H(1)^2+H(2)^2+H(3)^2);

% calculating semimajor a
a=-1*mu/(2*E);
% calculating essentricty ertyu
e=sqrt(1-(h^2/(mu*a)));

% calculating semiminor b
b=a*sqrt(1-e^2);
% calculating perigee & apogee rb ljb
rp=a*(1-e);
ra=a*(1+e);
% calculating period
period=2*pi*sqrt(a^3/mu);
g=acosd((1/e)*(((a*(1-e^2))/r)-(1)));

phi=acosd((1+e*cosd(g))/(sqrt(1 +2*e*cosd(g) +e^2 )));

 % % % % % % % % % % % % % % % 
 % % % % assignment two  % % % 
 % % % % % % % % % % % % % % % 

% calculating e vector ev
ev=(1/mu)*((v^2-(mu/r)).*R-(dot(R,V).*V));

% defining k and i unit vectors K I
K=[0 0 1];
I=[1 0 0];

% definig k.H KdH
KdH=dot(K,H);

% calculating inclination i
i=acosd(KdH/h);

% calculatin line of nodes N and N magnitud n
N=cross(K,H);
n=mag(N);
% calculating right ascensioin point omega
omega = acosd(dot(I,N)/n);
if N(2)<0
 omega =360-omega;
end

% calculating argument of perigee somega 
somega=acosd(dot(N,ev)/(n*e));
if ev(3)<0
    somega=360-somega;
end

% calculatin true anamoly g
g=acosd(dot(ev,R)/(e*r));
o=dot(V,R);
if o<0
    g=360-g;
end
% v
% r
% H

% E
% rp
% ra
% period 
% g
% b
a=44522.1
e=0;
i=36  ;
somega=0  ;       %+150  ;
omega=20  ;
period=2*pi*sqrt(a^3/mu);



somega=somega-90
omega=omega+90
i=180-i
% period=24*60*60
c=e*a/2;

axis([-a a -b b])
axis equal
rad = 6378.136;
subplot(2,1,1)       
axes(axes1)
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
camva(5);
  
%--------------------------------------------------

% plot equtorial plan -------------------------------
pointa1 = [-2*a,-2*b,0];
pointa2 = [-2*a,2*b,0];
pointb1 = [2*a,2*b,0];
pointb2 = [2*a,-2*b,0];
pointbs=[pointa1' pointa2' pointb1' pointb2']; % using the data given in the question
fill3(pointbs(1,:),pointbs(2,:),pointbs(3,:),'r')
grid on
alpha(0.5)
% ---------------------------------------------------


% % plot earth --------------------------------------6

earth=plotearth(0,0,0,rad);
direction = [0 0 1];
rotate(earth,direction,0)

axis equal
axis([-a a -b b -r r])
% % -------------------------------------------------
hold on
 
%  copied%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ni=150; %number of iterations
n=1; %nuber of roounds around the earth;
tr=100; %time ratio
dtime=period/ni;
M(1)=0;
E(1)=0;
pe2day=24*60*60/period;
for ii =1:ni *pe2day*2
%     time=time+dtime
    nn=sqrt(mu/a^3);
    M(ii+1)= dtime*nn+M(ii);
%     syms EE 
%     E(ii+1)= solve(EE-e*sin(EE) == M(ii+1),EE);
fun = @(EE) EE-e*sin(EE)-M(ii+1);  
x0 = E(ii); % initial point
E(ii+1) = fzero(fun,x0);
    nu(ii+1)=2*atan(tan(E(ii+1)/2)*sqrt((1+e)/(1-e)));
    if nu(ii+1) < 0
       nu(ii+1)=nu(ii+1)+ 2*pi;
    end
end


rad = 6378.136;
%  a = (rp + ra + 2*rad)/2;
      c = a - rp - rad;
%       e = c/a;
      p = a*(1 - e^2);
      n=1; %number of rotaions
      th = linspace(0,2*n*pi,70);
      th =nu;
      dthe=2*n*pi/70;
      r = p./(1 + e*cos(th));
      xx = r.*cos(th);
      yy = r.*sin(th);
      deg2rad=pi/180;
      omega = omega*deg2rad;
      i = i*deg2rad; 
      somega = somega*deg2rad; 
      axis tight
      
      % Coordinate Transformations
%       ZZ = [ cos(omega)   sin(omega)      0;
%              -sin(omega)  cos(omega)      0;
%                 0           0             1];
   ZZ = [ cos(i)   0      -sin(i);
           0           1             0;
          sin(i)    0     cos(i) ];
           
            
%       XX = [1          0             0;
%             0        cos(i)      sin(i);
%             0       -sin(i)     cos(i)];
      
      XX = [cos(somega)      sin(somega)     0;
            -sin(somega)     cos(somega)     0;
            0       0              1];
        
      ZZ2 = [ cos(omega)  sin(omega)    0;
             -sin(omega)  cos(omega)    0;
                 0           0            1];
             
      vec = XX*[xx;yy;zeros(1,length(xx))];
      vec=ZZ*vec;
      vec=ZZ2*vec;
      vec=fliplr(vec);


hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% actual plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
axis tight
% plot the orbit-----------------------------------------------------------
axis tight      
h6 = plot3(vec(1,:),vec(2,:),vec(3,:),'r');
    axis tight   

% -------------------------------------------------------------------------


% 2d image ****************************************************
subplot(2,1,2)       
axes(axes2)
axis tight
load topo
image([0 360],[-90 90], flip(topo), 'CDataMapping', 'scaled')
colormap(topomap1)

axis equal                                % set axis units to be the same size

ax = gca;                                 % get current axis
ax.XLim = [0 360];                        % set x limits
ax.YLim = [-90 90];                       % set y limits
ax.XTick = [0 60 120 180 240 300 360];    % define x ticks
ax.YTick = [-90 -60 -30 0 30 60 90];      % define y ticks
hold on
plot([0 360],[0 0],'r')
% ---------------------------------------------------------------------------

% motion plot &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subplot(2,1,1)    
% fake values to enter the loop with---------------------------------- 
   
h5 = plot3(0,0,0,'r*');
h6 = plot3(0,0,1,'r*');
 h7=plot3([0 ],[0 ],[0 3],'y*')
% ---------------------------------------------------------------------


% delit=360*period/day;
% % deli=-delit/length (vec);
deli=dtime/240
dsat=360/length (vec);
for t=1:length (vec)
        axes(axes1)
%         
         
          axis([-a a -b b -a a])
          x1(t)=vec(1,t);
          y1(t)=vec(2,t);
          z1(t)=vec(3,t);

        % drawing the satellite---------------------------------------
         delete(h5)
         h5 = plot3(x1(t),y1(t),z1(t),'r*');
% set(h5,'xdata',x1(t),'ydata',y1(t),'zdata',z1(t));
        % -------------------------------------------------------------

        % drawing line from center to satellite0----------------------
        delete(h6)
        h6 =  plot3([0 x1(t)],[0 y1(t)],[0 z1(t)],'color','yellow');
%           set(h6,'xdata',[0 x1(t)],'ydata',[0 y1(t)],'zdata',[0 z1(t)]);
        % --------------------------------------------------------------

        % drawing 3d ground track---------------------------------------
          R1=sqrt(x1(t)^2+y1(t)^2+z1(t)^2);
          x2(t)=(x1(t)/R1*rad)+sign(x1(t))*100;
          y2(t)=(y1(t)/R1*rad)+sign(y1(t))*100;
          z2(t)=(z1(t)/R1*rad)+sign(z1(t))*100;
           delete(h7)
          h7=plot3([0 x2],[0 y2],[0 z2],'y*');
        % --------------------------------------------------------------

        
        % rotation of camera with ground track--------------------------
%           view((th(t)/pi*180)+180,20);
          [azi(t),ele(t),r(t)] = cart2sph(x1(t),y1(t),z1(t));
           eli(t)=ele(t)*180/pi;
           aze(t)=azi(t)*180/pi;
           view(-(-aze(t)+90)+180+20,20)
           
        %   --------------------------------------------------------------
axis off
% camva(1)
% view(0,90)

        % 2dground tracks--------------------------------------------
       subplot(2,1,2)       
        axes(axes2)
        hold on
        axis([0 360 -90 90])
        hold on
        aaa=(aze);
while aaa(t)-(t*deli)<0
      aaa(t)=aaa(t)+360;
end

while aaa(t)-(t*deli)> 360
      aaa(t)=aaa(t)-360;
end

        u(t)=aaa(t)-(t*deli);
        GT=  plot(aaa(t)-(t*deli),-eli(t),'y*');     
        drawnow 
% 

% rotating earth model---------------------------------------------
% % 
% deli=dthe*180/pi;
 ZZ3 = [ cosd(deli)  -sind(deli)    0;
        sind(deli)  cosd(deli)    0;
         0           0           1];
             
 gt=[ x2;y2;z2];
 ngt=ZZ3*gt;
 x2=ngt(1,:);
 y2=ngt(2,:);
 z2=ngt(3,:);

direction = [0 0 1];
% rotate(h7,direction,deli)
rotate(earth,[0 0 1],deli,[0 0 0]) 
pause(dtime/100000)
% pause(.00003)
end
       subplot(2,1,1)       
       axes(axes1)
%        h7=plot3([0 x2],[0 y2],[0 z2],'y*');
       rotate3d on
% nu
% size(nu)

