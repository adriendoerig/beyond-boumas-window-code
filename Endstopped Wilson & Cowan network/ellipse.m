function [x,y] = ellipse(n,a,b,h,k,theta)

% ELLIPSE create an ellipse in the plane.

% The equation of an ellipse is:

%

% ( (x-h)/a ).^2 + ( (y-k)/b ).^2 = 1

%

% Where:

%

% x = the x-coordinate

% y = the y-coordinate

% h = the displacement of the center of the ellipse along the x-axis

% k = the displacement of the center of the ellipse along the y-axis

% 2*a = the length of the ellipse

% 2*b = the height of the ellipse

%

% This function solves for y as a function of x, a, b, h, and k.

%

% Usage:

%

% [x,y] = ellipse(n,a,b,h,k,theta);

%

% n = the number of elements in the vectors x and y

% E = the rotation matrix

%

% e.g.

% n = 100;

% a = 10;

% b = 5;

% h = -2;

% k = 3;

% Theta = pi/4;

% [x,y] = ellipse(n,a,b,h,k,Theta);

% figure('Color',[1 1 1]);

% plot(x,y,'k-',h,k,'go','LineWidth',2);

%

% Created by Aaron Clarke, York University, Elder Lab, October 12, 2003

if a<=0 || b<=0

    x = h; y = k;

    return

end

x1 = linspace((h-a),(h+a),n);

y1 = real(k + sqrt(b^2*(1-(x1.^2 - 2*x1*h + h^2)/(a^2))));

x2 = fliplr(x1); %x1; % Remove comment for line through ellipse.

y2 = real(k - sqrt(b^2*(1-(x2.^2 - 2*x2*h + h^2)/(a^2))));

y = [y1 y2];

x = [x1 x2];

if theta~=0

    theta = -theta;

    x = x-h;

    y = y-k;       

    

    u = cos(theta)*x + sin(theta)*y;

    v = -sin(theta)*x + cos(theta)*y;

   

    x = u+h;

    y = v+k;

end

   

 

return