function S=distance_contextual_elements(maxx,maxy,distance,nr)
%      |||||
% |||||||||||||||
% |||||||||||||||
%      |||||
% vertical distance of contextual elements (nr pixels)
% nr = number of contextual elements
% max distance=19



h=29; %hoogte vernier element
h2=19; %hoogte contextual element
d=2; %halve afstand tussen twee vernier elementen

S = zeros (maxx,maxy);
cx=round(maxx/2);
cy=round(maxy/2);

r=cx-(nr-1)/2*10:10:cx+(nr-1)/2*10;

%grating elements
S(cx-120:10:cx+120,[cy-h-d:cy-d,cy+d:cy+h+d]) = 1;

%contextual elements
v1=cy-d-h-distance-h2:cy-d-h-distance;
v2=cy+d+h+distance:cy+d+h+distance+h2;

S(r,[v1,v2]) = 1;