function[e]= plotearth(x,y,z,rad)

          [X,Y,Z] = sphere(50);
      load topo
      topo = [topo(:,181:360) topo(:,1:180)];
      mat.dull.AmbientStrength = 0.4;
      mat.dull.DiffuseStrength = .6;
      mat.dull.SpecularColorReflectance = .5;
      mat.dull.backfacelighting = 'reverselit';
      mat.dull.SpecularExponent = 20;
      mat.dull.SpecularStrength = .8;
      e=surface((rad*X)+x,(rad*Y)+y,(rad*Z)+z, ...
         mat.dull, ...
         'FaceColor','texturemap',...
         'EdgeColor','none',...
         'FaceLighting','phong',...
         'Cdata',topo,'tag','earth');   
      colormap(topomap1)
      light('position',rad*[10 10 10]);
      %light('position',rad*[-10 -10 -10], 'color', [.6 .2 .2]);
      %set(gcf,'renderer','opengl');
           axis vis3d