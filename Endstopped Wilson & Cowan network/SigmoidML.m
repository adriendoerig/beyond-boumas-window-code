function s = SigmoidML(parms,x)

% SIGMOIDML Sigmoidal function for maximum likelihood function fitting.
% 
% Usage: s = SigmoidML(parms,x);
% 
% e.g.
% 
% parms = [50 0.5 2 4];
% x = 0:1:100;
% s = SigmoidML(parms,x);
% 
% figure;
% plot(x,s,'-r');
% 
% Created by Aaron Clarke, NDSU, Visual Dynamics lab, July 23, 2008

if sum(parms)==0,
   s=x;
   return
end

m = parms(1);
slope = parms(2);
a1 = parms(3);
a2 = parms(4);

%if h<0,
%   s=zeros(size(x));
%   return
%end

s = (a2-a1)./(exp(-(x-m)*slope)+1) + a1;
