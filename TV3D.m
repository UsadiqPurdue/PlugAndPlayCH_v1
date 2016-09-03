function [TV]=TV3D(x)
   
   [g1,g2,g3]=gradient(x);
   
   g=[g1(:) g2(:) g3(:)];

   TV=norm(g,1);
   
end