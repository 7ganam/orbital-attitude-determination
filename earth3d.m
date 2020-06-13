clf
[x,y,z] = sphere(50);          % create a sphere
s = surface(x,y,z); 
load topo    % plot spherical surfacec
s.CData = topo;                % set color data to topographic data
s.FaceColor = 'texturemap';    % use texture mapping
s.EdgeColor = 'none';          % remove edges
s.FaceLighting = 'gouraud';    % preferred lighting for curved surfaces
s.SpecularStrength = 0.4;      % change the strength of the reflected light

light('Position',[-1 0 1])     % add a light

axis square off                % set axis to square and remove axis
view([-30,30])                 % set the viewing angle
% axis on
