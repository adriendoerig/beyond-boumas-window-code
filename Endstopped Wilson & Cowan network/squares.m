function S=squares(size)

d=2;
% vernier
   S = zeros (300,140);
   S(150-d,39:68) = 1;
   S(150+d,72:101) = 1;
   
%squares
   S(140,39:101) = 1;
   S(160,39:101) = 1;
   S(140-size:140,39) = 1;
   S(140-size:140,101) = 1;
   S(160:160+size,39) = 1;
   S(160:160+size,101) = 1;
   S(140-size,39:101) = 1;
   S(160+size,39:101) = 1;
   

   