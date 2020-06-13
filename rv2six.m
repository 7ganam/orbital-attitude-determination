function[i,omega,a,h,phi,somega,g,e,rp,ra,b,E]= rv2six(R,V) 
v=sqrt(V(1)^2+V(2)^2+V(3)^2);
r=sqrt(R(1)^2+R(2)^2+R(3)^2);
 global mu;
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