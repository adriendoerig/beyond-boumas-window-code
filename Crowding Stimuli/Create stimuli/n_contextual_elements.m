function S=n_contextual_elements(maxx,maxy,nr)
%      |||||
% |||||||||||||||
% |||||||||||||||
%      |||||


n=(nr-1)/2;

h=29; %hoogte vernier element
h2=19; %hoogte contextual element
d=2; %halve afstand tussen twee vernier elementen
d2=11; %afstand grating elements and contextual elements

S = zeros (maxx,maxy);
cx=round(maxx/2);
cy=round(maxy/2);


%grating elements
S(cx-120:10:cx+120,[cy-h-d:cy-d,cy+d:cy+h+d]) = 1;

%contextual elements
v1=cy-d-h-h2-d2:cy-d-h-d2;
v2=cy+d+h+d2:cy+d+h+h2+d2;

S([cx-n*10:10:cx+n*10],[v1,v2]) = 1;