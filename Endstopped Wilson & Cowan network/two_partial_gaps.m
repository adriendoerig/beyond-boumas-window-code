function S=two_partial_gaps(pos,L);

d=29-L;

S = zeros (300,140);
S(30:10:270,[39:68,72:101]) = 1;
if d>=0,
    S([150-pos*10,150+pos*10],[39:39+d,101-d:101]) = 0;
elseif d<0,    
    S([150-pos*10,150+pos*10],[39+d:68,72:101-d]) = 1;
end;    