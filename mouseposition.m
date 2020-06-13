function[x,y,z,rad]= plotearth()


load coast
         ncst = ncst *pi/180;
         all =zeros(length(ncst),3);
         for j = 1:length(ncst)
            theta = ncst(j,1);
            phi = ncst(j,2);
            all(j,:) = [cos(theta)*cos(phi),...
                        sin(theta)*cos(phi),...
                        -sin(phi)];
         end
         rad=152
        plot3(rad*all(:,1),rad*all(:,2),-rad*all(:,3));
%          set(data.earthhandle(2),'color',[0 .9 0]);
          [X,Y,Z] = sphere(50);
      load topo
      topo = [topo(:,181:360) topo(:,1:180)];
      mat.dull.AmbientStrength = 0.4;
      mat.dull.DiffuseStrength = .6;
      mat.dull.SpecularColorReflectance = .5;
      mat.dull.backfacelighting = 'reverselit';
      mat.dull.SpecularExponent = 20;
      mat.dull.SpecularStrength = .8;
      surface(rad*X,rad*Y,rad*Z, ...
         mat.dull, ...
         'FaceColor','texturemap',...
         'EdgeColor','none',...
         'FaceLighting','phong',...
         'Cdata',topo,'tag','earth');   
      colormap(topomap1)
      light('position',rad*[10 10 10]);
      %light('position',rad*[-10 -10 -10], 'color', [.6 .2 .2]);
      %set(gcf,'renderer','opengl');
      data.simple = 1;