function S=Xsquares(size)

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
   
%crosses
   x = [140-size, 140];
   y = [39, 101];
   l = round(interp1(x,y, 140-size:140));
   for i = 1:size
       S(140-size-1+i, l(i)-1:l(i):l(i)+1) = 1;
   end
   
   x = [140-size, 140];
   y = [101, 39];
   l = round(interp1(x,y, 140-size:140));
   for i = 1:size
       S(140-size-1+i, l(i)-1:l(i):l(i)+1) = 1;
   end
   
   x = [160, 160+size];
   y = [101, 39];
   l = round(interp1(x,y, 160:160+size));
   for i = 1:size
       S(160-1+i, l(i)-1:l(i):l(i)+1) = 1;
   end
   
   x = [160, 160+size];
   y = [39, 101];
   l = round(interp1(x,y, 160:160+size));
   for i = 1:size
       S(160-1+i, l(i)-1:l(i):l(i)+1) = 1;
   end
   