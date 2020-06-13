function varargout = orbital_project(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @orbital_project_OpeningFcn, ...
                   'gui_OutputFcn',  @orbital_project_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before orbital_project is made visible.
function orbital_project_OpeningFcn(hObject, eventdata, handles, varargin)
zoo=.7;
% set(handles.numb,'string',1)
set(handles.x,'string',8228)
set(handles.y,'string',389)
set(handles.z,'string',6888)
set(handles.vx,'string',-.7)
set(handles.vy,'string',6.6)
set(handles.vz,'string',-.6)
rotate3d off
 figure(1)       
 axes1=gca;
 f1=gcf;
 set(f1, 'Position',[ 73   136   275   241])
 axis equal
 axis off
 axis vis3d
 camva(.7)

 figure(2)     
 axes2=gca;
 f2=gcf
 set(f2, 'Position',[368   142   275   234])
 camva(10)
 figure(3)       
 axes3=gca;
 f3=gcf;
 set(f3, 'Position',[680   136   276   237])
 
 figure(4)     
 axes4=gca;
 f4=gcf
 set(f4, 'Position',[976   134   275   241])
 

handles.output = hObject;
guidata(hObject, handles);
function varargout = orbital_project_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
function x_Callback(hObject, eventdata, handles)
function x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in run.
%%rotating Earth code RV
function run_Callback(hObject, eventdata, handles)
 rotate3d off

 global mu;
 mu=3.9868*10^5;
 blah=get(handles.came,'selectedobject');
 came=get(blah,'string');
 came=strcmp(came,'rotation');
 
 bl=get(handles.sit,'selectedobject');
 sit=get(bl,'string');
 
 
 figure(1)      
 axes1=gca;
 axis equal
 axis off
 axis vis3d
 camva(.7)
 
 figure(2)     
 axes2=gca;

 camva(8)
 
 figure(3)       
 axes3=gca;

 figure(4)     
 axes4=gca;

 
cla(axes3,'reset');
cla(axes4,'reset');
%% getting R and V vectors from the GUI
x=str2double(get(handles.x,'string'));
y=str2double(get(handles.y,'string'));
z=str2double(get(handles.z,'string'));
vx=str2double(get(handles.vx,'string'));
vy=str2double(get(handles.vy,'string'));
vz=str2double(get(handles.vz,'string'));

V=[vx vy vz];
R= [x y z];
% % % %calculating 6oe 
[v0 r0 h0 a0 b0 e0 E0 rp g0 somega0 omega0 i0]=rotatingEart_RV(R,V,axes3,axes4,came)
%% setting output to the GUI
set(handles.v,'string',v0)
set(handles.r,'string',r0)
set(handles.h,'string',h0)
set(handles.a,'string',a0)
set(handles.b,'string',b0)
set(handles.e,'string',e0)
set(handles.E,'string',E0)
set(handles.rp,'string',rp0)
set(handles.g,'string',g0)
set(handles.somega,'string',somega0)
set(handles.omega,'string',omega0)
set(handles.i,'string',i0)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 




% non rotating RV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nrun_Callback(hObject, eventdata, handles)
%
blah=get(handles.came,'selectedobject');
 came=get(blah,'string');
 came=strcmp(came,'rotation');
rotate3d off
  figure(1)       
 axes1=gca;

 axis equal
 axis off
 axis vis3d
 camva(.7)
 
 figure(2)     
 axes2=gca;
 
 camva(10)
 figure(3)       

 
 figure(4)     

 
cla(axes1,'reset');
cla(axes2,'reset');





rotate3d off


x=str2double(get(handles.x,'string'));
y=str2double(get(handles.y,'string'));
z=str2double(get(handles.z,'string'));
vx=str2double(get(handles.vx,'string'));
vy=str2double(get(handles.vy,'string'));
vz=str2double(get(handles.vz,'string'));

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

[v0 r0 h0 a0 b0 e0 E0 rp0 g0 somega0 omega0 i0]=non_rotatingEart_RV(R,V,axes1,axes2,came)
%% setting output to the GUI
set(handles.v,'string',v0)
set(handles.r,'string',r0)
set(handles.h,'string',h0)
set(handles.a,'string',a0)
set(handles.b,'string',b0)
set(handles.e,'string',e0)
set(handles.E,'string',E0)
set(handles.rp,'string',rp0)
set(handles.g,'string',g0)
set(handles.somega,'string',somega0)
set(handles.omega,'string',omega0)
set(handles.i,'string',i0)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 




% % % % % % % % % % % % % % % % nonrotating coe
% --- Executes on button press in NONrot.
function NONrot_Callback(hObject, eventdata, handles)
blah=get(handles.came,'selectedobject');
 came=get(blah,'string');
 came=strcmp(came,'rotation');
rotate3d off
  figure(1)       
 axes1=gca;

 axis equal
 axis off
 axis vis3d
 camva(.7)
 
 figure(2)     
 axes2=gca;
 
 camva(10)
 figure(3)       

 
 figure(4)     

 
cla(axes1,'reset');
cla(axes2,'reset');





rotate3d off



% assignement main code
% input V & R 
mu=3.9868*10^5;



% x=str2double(get(handles.x,'string'));
a=str2double(get(handles.a,'string'));
e=str2double(get(handles.e,'string'));
i=str2double(get(handles.i,'string'));
somega=str2double(get(handles.somega,'string'));
omega=str2double(get(handles.somega,'string'));
g=str2double(get(handles.g,'string'));

global mu;
h=sqrt(a*(1-e^2)*mu);

coe(1)=h;
coe(2)=e;
coe(3)=omega;
coe(4)=i;
coe(5)=somega;
coe(6)=g;
 [R, V] = sv_from_coe(coe,mu);

% find the magnitude of R and V
v=sqrt(V(1)^2+V(2)^2+V(3)^2);
r=sqrt(R(1)^2+R(2)^2+R(3)^2);

[v0 r0 h0 a0 b0 e0 E0 rp0 g0 somega0 omega0 i0]=non_rotatingEart_RV(R,V,axes1,axes2,came)
%% setting output to the GUI
set(handles.v,'string',v0)
set(handles.r,'string',r0)
set(handles.h,'string',h0)
set(handles.a,'string',a0)
set(handles.b,'string',b0)
set(handles.e,'string',e0)
set(handles.E,'string',E0)
set(handles.rp,'string',rp0)
set(handles.g,'string',g0)
set(handles.somega,'string',somega0)
set(handles.omega,'string',omega0)
set(handles.i,'string',i0)



% --- Executes on button press in ROT.
function ROT_Callback(hObject, eventdata, handles)

mu=3.9868*10^5;
rotate3d off
 zoo=.7;
 blah=get(handles.came,'selectedobject');
 came=get(blah,'string')
 came=strcmp(came,'rotation')
 
 bl=get(handles.sit,'selectedobject');
 sit=get(bl,'string')
 
 
 figure(1)       
 axes1=gca;
 f1=gcf;
%  set(f1, 'Position',[ 55   213   307   275])
 axis equal
 axis off
 axis vis3d
 camva(.7)
 
 figure(2)     
 axes2=gca;
 f2=gcf
%  set(f2, 'Position',[ 379   213   309   275])
 camva(8)
 figure(3)       
 axes3=gca;
 f3=gcf;
%  set(f3, 'Position',[704   214   307   273])
 
 figure(4)     
 axes4=gca;
 f4=gcf
%  set(f4, 'Position',[1029         223         304         263])
cla(axes3,'reset');
cla(axes4,'reset');
rotate3d off
a=str2double(get(handles.a,'string'));
e=str2double(get(handles.e,'string'));
i=str2double(get(handles.i,'string'));
somega=str2double(get(handles.somega,'string'));
omega=str2double(get(handles.somega,'string'));
vz=str2double(get(handles.vz,'string'));
% hObject    handle to ROT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)% hObject    handle to ROT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% a=44522.1
% e=0;
% i=36  ;
% somega=0  ;       %+150  ;
% omega=20  ;
% period=2*pi*sqrt(a^3/mu);

% buffering to sync with stk
period=2*pi*sqrt(a^3/mu);
somega=somega-90;
omega=omega+90;
i=180-i;

% period=24*60*60






c=e*a/2;
%----------------------------

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
% axis([-a a -b b])
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



hold on
 
%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate time dependant matrix-------------------------------------------
ni=200; %number of iterations
n=1; %nuber of roounds around the earth;
tr=100; %time ratio
dtime=period/ni;
M(1)=0;
E(1)=0;
pe2day=24*60*60/period;


for ii =1:ni%*pe2day
  nn=sqrt(mu/a^3);
    M(ii+1)= dtime*nn+M(ii);
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
%       th = linspace(0,2*n*pi,70);
      th =nu;
      dthe=2*n*pi/70;
      r = p./(1 + e*cos(th));
      xx = r.*cos(th);
      yy = r.*sin(th);
      deg2rad=pi/180;
      omega = omega*deg2rad;
      i = i*deg2rad; 
      somega = somega*deg2rad; 
      
      
      % Coordinate Transformations

      ZZ = [ cos(i)   0      -sin(i);
              0           1             0;
             sin(i)    0     cos(i) ];
           
        
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

% plot the orbit-----------------------------------------------------------
     h6 = plot3(vec(1,:),vec(2,:),vec(3,:),'r');
% -------------------------------------------------------------------------


% 2d image ----------------------------------------------------------------
figure(4)       
axes(axes4)
% axis tight
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
plot([0 360],[0 0],'b')
% ---------------------------------------------------------------------------


% motion plot &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
figure(3)    
 view(30 ,20)
% fake values to enter the loop with----------------------------------  
h5 = plot3(0,0,0,'r*');
h6 = plot3(0,0,1,'r*');
h7=plot3([0 ],[0 ],[0 3],'y*');
% ---------------------------------------------------------------------

camva(8)
deli=dtime/240;
figure(3)
axis equal
axis off
axis vis3d
camlookat(findobj('tag','earth'));
camva(.7)
 view(30 ,20)

%% the main loop
for t=1:length (vec)
        axes(axes3)         
          x1(t)=vec(1,t);
          y1(t)=vec(2,t);
          z1(t)=vec(3,t);
     
        % drawing the satellite---------------------------------------
         delete(h5)
         h5 = plot3(x1(t),y1(t),z1(t),'r*');
% set(h5,'xdata',x1(t),'ydata',y1(t),'zdata',z1(t));
        % -------------------------------------------------------------
% axis ([-a a -a a])
        % drawing line from center to satellite0----------------------
        delete(h6)
        h6 =  plot3([0 x1(t)],[0 y1(t)],[0 z1(t)],'color','yellow');
%           set(h6,'xdata',[0 x1(t)],'ydata',[0 y1(t)],'zdata',[0 z1(t)]);
        % --------------------------------------------------------------
      axis vis3d
        % drawing 3d ground track---------------------------------------
          R1=sqrt(x1(t)^2+y1(t)^2+z1(t)^2);
          x2(t)=(x1(t)/R1*rad)+sign(x1(t))*100;
          y2(t)=(y1(t)/R1*rad)+sign(y1(t))*100;
          z2(t)=(z1(t)/R1*rad)+sign(z1(t))*100;
           delete(h7)
          h7=plot3([0 x2],[0 y2],[0 z2],'o','MarkerEdgeColor','yellow','MarkerFaceColor','yellow','MarkerSize', 3);
        % --------------------------------------------------------------
% plot (long(t),lat(t),'o', 'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize', 3);
        
        % rotation of camera with ground track--------------------------
%           view((th(t)/pi*180)+180,20);y

          [azi(t),ele(t),r(t)] = cart2sph(x1(t),y1(t),z1(t));
           eli(t)=ele(t)*180/pi;
      
          aze(t)=azi(t)*180/pi;
          if came ==1
          view(-(-aze(t)+90)+180+20,20)
          else
              view(30 ,20)
          end
%    --------------------------------------------------------------


        % 2dground tracks--------------------------------------------
       figure(4)
        axes(axes4)
        hold on
        axis([0 360 -90 90])
        hold on
        %make sure the angels within the 2d image borders-----------------
        aaa=(aze);
            while aaa(t)-(t*deli)<0
                  aaa(t)=aaa(t)+360;
            end

            while aaa(t)-(t*deli)> 360
                  aaa(t)=aaa(t)-360;
            end
         %-----------------------------------------------------------------
                   long(t)=aaa(t)-(t*deli);
                   lat(t)=-eli(t);
%plot 2d ground track-----------------------------------------------------
        GT=  plot(aaa(t)-(t*deli),-eli(t),'o','MarkerEdgeColor', 'yellow','MarkerFaceColor','yellow','MarkerSize', 3);
%     
%     plot (long(t),lat(t),'o', 'MarkerEdgeColor', 'k','MarkerFaceColor','r','MarkerSize', 3);
%     if (t~=1 && abs(long(t-1)-long(t))<100)
%         line([long(t-1) long(t)],[lat(t-1) lat(t)],'Color', 'red', 'LineWidth', 2);
%     end




% ------------------------------------------------------------------------

% rotating earth model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  find new values for 3d ground track--------------------------------
ZZ3 = [ cosd(deli)  -sind(deli)    0;
        sind(deli)  cosd(deli)    0;
         0           0           1];
 gt=[ x2;y2;z2];
 ngt=ZZ3*gt;
 x2=ngt(1,:);
 y2=ngt(2,:);
 z2=ngt(3,:);
%  -----------------------------------------------------------------------

rotate(earth,[0 0 1],deli,[0 0 0]) 
pause(dtime/100000)
end

% --- Executes during object creation, after setting all properties.
function numb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edit67_Callback(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit67 as text
%        str2double(get(hObject,'String')) returns contents of edit67 as a double


% --- Executes during object creation, after setting all properties.
function edit67_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R12_Callback(hObject, eventdata, handles)
% hObject    handle to R12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R12 as text
%        str2double(get(hObject,'String')) returns contents of R12 as a double


% --- Executes during object creation, after setting all properties.
function R12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R13_Callback(hObject, eventdata, handles)
% hObject    handle to R13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R13 as text
%        str2double(get(hObject,'String')) returns contents of R13 as a double


% --- Executes during object creation, after setting all properties.
function R13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R21_Callback(hObject, eventdata, handles)
% hObject    handle to R21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R21 as text
%        str2double(get(hObject,'String')) returns contents of R21 as a double


% --- Executes during object creation, after setting all properties.
function R21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R22_Callback(hObject, eventdata, handles)
% hObject    handle to R22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R22 as text
%        str2double(get(hObject,'String')) returns contents of R22 as a double


% --- Executes during object creation, after setting all properties.
function R22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R23_Callback(hObject, eventdata, handles)
% hObject    handle to R23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R23 as text
%        str2double(get(hObject,'String')) returns contents of R23 as a double


% --- Executes during object creation, after setting all properties.
function R23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R31_Callback(hObject, eventdata, handles)
% hObject    handle to R31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R31 as text
%        str2double(get(hObject,'String')) returns contents of R31 as a double


% --- Executes during object creation, after setting all properties.
function R31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R32_Callback(hObject, eventdata, handles)
% hObject    handle to R32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R32 as text
%        str2double(get(hObject,'String')) returns contents of R32 as a double


% --- Executes during object creation, after setting all properties.
function R32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R33_Callback(hObject, eventdata, handles)
% hObject    handle to R33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R33 as text
%        str2double(get(hObject,'String')) returns contents of R33 as a double


% --- Executes during object creation, after setting all properties.
function R33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
