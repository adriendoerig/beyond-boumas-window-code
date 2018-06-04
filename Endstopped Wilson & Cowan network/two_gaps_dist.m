function S=two_gaps_dist(dist);

S = zeros (300,140);
S(130:10:170,[39:68,72:101]) = 1;
S([40-dist:10:130-dist,170+dist:10:260+dist],[39:68,72:101]) = 1;