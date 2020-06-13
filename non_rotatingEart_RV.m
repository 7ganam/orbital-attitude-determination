function [v0 r0 h0 a0 b0 e0 E0 rp0 g0 somega0 omega0 i0]=non_rotatingEart_RV(R,V,axes1,axes2,came)
E0=0;
global mu;
%%calculating 6OE
% find the magnitude of R and V
v=sqrt(V(1)^2+V(2)^2+V(3)^2);
r=sqrt(R(1)^2+R(2)^2+R(3)^2);
rad2degr=180/pi;
% 
% [i,omega,a,h,phi,somega,g,e,rp,ra,b,E]= rv2six(R,V) 
[h e omega i somega g a rp ra b E]= coe_from_sv(R,V,mu);
omega=omega*rad2degr;
i=i*rad2degr;
somega =somega *rad2degr;
g=g*rad2degr;
%%
%% setting output to the GUI
v0=v;
r0=r;
h0=h;
a0=a;
b0=b;
e0=e;
E0=E;
rp0=rp;
g0=g;
somega0=somega;
omega0=omega;
i0=i;
 

% buffering to sync with stk
period=2*pi*sqrt(a^3/mu);
somega=somega-90;
omega=omega+90;
i=180-i;
%----------------------------

earth=preplot(axes1,a,b);
hold on
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate time dependant matrix-------------------------------------------
ni=200; %number of iterations
n=1; %nuber of roounds around the earth;
dtime=period/ni;
g
sqrt( (1-e)/(1+e))* tand(g/2) 
E0= 2*atan(sqrt( (1-e)/(1+e))* tand(g/2) );
if E0<0
    E0=E0+2*pi;
end
M0=E0-e*sin(E0);
t0=M0/(2*pi)*period

fprintf('E)= %d',E0);
M(1)=M0;
E(1)=E0;

for ii =1:ni%*pe2day
    tf=t0+dtime*ii;
    if tf>period
        tf=tf-period;
    end
    
    nn=sqrt(mu/a^3);
    M(ii+1)= tf*nn;
    fun = @(EE) EE-e*sin(EE)-M(ii+1);  
    x0 = E(ii); % initial point
    E(ii+1) = fzero(fun,x0);
    nu(ii+1)=2*atan(tan(E(ii+1)/2)*sqrt((1+e)/(1-e)));
    if nu(ii+1) < 0
       nu(ii+1)=nu(ii+1)+ 2*pi;
    end
end

%------------------------------------------------------------------------


      rad = 6378.136;
      p = a*(1 - e^2);
      n=1; %number of rotaions
      th =nu;
      r = p./(1 + e*cos(th));
      xx = r.*cos(th);
      yy = r.*sin(th);
      deg2rad=pi/180;
      omega = omega*deg2rad;
      i = i*deg2rad; 
      somega = somega*deg2rad; 
      
j2=1.08263 * 10^(-3);
size(r)
for ii=1:length(xx)
    
    %%earth oblateness correction
omegad= -(3/2)*( (sqrt(mu)*j2*(r(ii)^2))  /  ((1-e)^2 * a^(7/2) )   )*cos(i);
somegad= -(3/2)*( (sqrt(mu)*j2*r(ii)^2)  /  ((1-e)^2 * a^(7/2) )   )*( (5/2)*sin(i)*sin(i)-2 );

omega_new=omega+omegad;
somega_new=somega+somegad;
      % Coordinate Transformations

      ZZ = [ cos(i)   0      -sin(i);
              0           1             0;
             sin(i)    0     cos(i) ];
           
        
      XX = [cos(somega_new)      sin(somega_new)     0;
            -sin(somega_new)     cos(somega_new)     0;
            0       0              1];
        
      ZZ2 = [ cos(omega_new)  sin(omega_new)    0;
             -sin(omega_new)  cos(omega_new)    0;
                 0           0            1];
             
       vec(:,ii) = XX*[xx(ii);yy(ii);0];
       vec(:,ii)=ZZ*  vec(:,ii);
       vec(:,ii)=ZZ2*  vec(:,ii);
       vec(:,ii)=fliplr(  vec(:,ii));
hold on
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% actual plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
% plot the orbit-----------------------------------------------------------
     
h6 = plot3(vec(1,:),vec(2,:),vec(3,:),'r');
%     axis tight   
% -------------------------------------------------------------------------


% 2d image ****************************************************
    axes(axes2)
%     axis tight
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

% fake values to enter the loop with---------------------------------- 
h5 = plot3(0,0,0,'r*');
h6 = plot3(0,0,1,'r*');
% ---------------------------------------------------------------------
camva(8)
figure(1)
axis equal
 axis off
 axis vis3d
 camva(.7)

 
 
 
for t=1:length (vec)
        
         axes(axes1)
          delete(h5)
          delete(h6)

          x1=vec(1,length(vec)-t+1);
          y1=vec(2,length(vec)-t+1);
          z1=vec(3,length(vec)-t+1);
        % drawing the satellite---------------------------------------
        h5 = plot3(x1,y1,z1,'r*');
        % -------------------------------------------------------------

        % drawing line from center to satellite0----------------------
        h6 =  plot3([0 x1],[0 y1],[0 z1],'color','yellow');
        % --------------------------------------------------------------

        % drawing 3d ground track---------------------------------------
          R1=sqrt(x1^2+y1^2+z1^2);
          x2=(x1/R1*rad)+sign(x1)*100;
          y2=(y1/R1*rad)+sign(y1)*100;
          z2=(z1/R1*rad)+sign(z1)*100;
          h7=plot3([0 x2],[0 y2],[0 z2],'y*');
        % --------------------------------------------------------------

        
        % rotation of camera with ground track--------------------------
%           view((th(t)/pi*180)+180,20);
          [azi(t),ele(t),r(t)] = cart2sph(x1,y1,z1);
           eli(t)=ele(t)*180/pi;
           aze(t)=azi(t)*180/pi;
          if came ==1
          view(-(-aze(t)+90)+180+20,20)
          else
              view(30 ,20)
          end
        %   --------------------------------------------------------------

        % 2dground tracks--------------------------------------------
     
        axes(axes2)
        hold on
        axis([0 360 -90 90])
        aaa=(aze);
        if aaa(t)<0
           aaa(t)=aaa(t)+360;
            
        end
        long(t)=aaa(t);
        lat(t)=-eli(t);
               
    plot (long(t),lat(t),'o', 'MarkerEdgeColor', 'yellow','MarkerFaceColor','yellow','MarkerSize', 3);

                  drawnow
end