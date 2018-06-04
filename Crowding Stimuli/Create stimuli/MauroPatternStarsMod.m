function [Img,Data,XCoords,YCoords] = MauroPatternStars(Cond,side)

% MAUROPATTERNSTARS Generate the stimuli from Mauro's "pattern_irregular"
% crowding experiment.
%
% Usage: [Img,Data,XCoords,YCoords] = MauroPatternStars(Condition);
%
% Conditions may be any integer from 1 to 14
%
% e.g.
%
% for Cond = 1:14
%     [Img,Data,XCoords,YCoords] = MauroPatternStars(Cond);
%     for iVernOffst = 1:size(Img,3)
%         figure;
%         imshow(Img(:,:,iVernOffst));
%         figure('Position',[19 619 614 296]);
%         plot(XCoords,YCoords,'k-');
%         axis equal;
%         axis([0 450 0 210]);
%         axis ij;
%         set(gca,'XTick',[],'YTick',[]);
%     end
% end
%
% Created by Mauro Manassi (Experiment_pattern_stars.m)
% Modified by Aaron Clarke, 27 October 2014, EPFL, LPSY Herzog Lab to only
% produce the stimuli and not run the experiment.

data = [139.5	1497	502	1993	98.5	1994	1993	1489.333333	1989.5;
    41	1290	269	1405	37.5	1367.5	1994.5	85	1020.5;
    96.5	1986.5	954.5	1576	857.6666667	1219	1619.5	1625	1439;
    212	1212.5	633.3333333	1937.5	910	1990	1968.5	1891.5	1918;
    42.5	1942	937.3333333	1984.5	737.3333333	1943	1997.5	1976	1972;
    128	1990	349.5	1995.5	517	1988	1873	1919	1945.5]; % We're not sure if data(:,10) is the data for stimulus 10
Data = zeros(2,length(data(1,:)));
Data(1,:) = sum(data,1)/length(data(:,1));
Data(2,:) = std(data,1)/sqrt(length(data(:,1))); 
clarkeFactor = 1;% just to check if Aaron's code works only because of the specific stim size he chose. Set to one to cancel effects.

% Stimulus conversion factors
SecPerPix = 275;                    % arc seconds per pixel
SecPerDeg = 3600;                   % seconds of arc per degree of visual angle
PixPerSec = 1/SecPerPix;            % Pixels per second of arc
PixPerMin = 60/SecPerPix;           % Pixels per minute of arc
PixPerDeg = SecPerDeg/SecPerPix;    % Pixels per degree of visual angle

Offst = 1.5;                        % Vernier offset (pix)
ECC_FIX = 9;                        % [deg] %from the center
VERNIER_LEN = 2400;                 % [arcsec], one segment
VERNIER_GAP = 240;                  % [arcsec]
POINTER_GAP = 18000;                % [arcsec] for the pointers
ECC_VER = 0;
FLANKER_LEN = 7200;                 % [arcsec]
FIX_DIA = 8;                        % [arcmin]
sdLine = 0.5;                       % standard error of the gaussian,
                                    % corresponds to line width
freq = 10; % spacing between flankers

INT_FLANKER = 3600;                 % [arcsec]
IRR_R = 1.4*(FLANKER_LEN/2);

%-----------

% Image dimensions
height = 210;
width = 450;
ImSz = [height width];
Img = zeros(height,width,1);

Xcenter = 0.5*width;   % center X coordinates
Ycenter = 0.5*height;  % center Y coordinates

SPACING = 2.7;                  % [deg]  %from the center for the vernier
Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing

SPACING = 2.2;% [deg]  %from the center for the vernier
ECC_VER1 = ECC_VER+SPACING; % [deg]  %from the center for square 1
ECC_VER2 = ECC_VER-SPACING;% [deg]  %from the center for square 2
ECC_VER3 = ECC_VER+2*SPACING;% [deg]  %from the center for square 3
ECC_VER4 = ECC_VER-2*SPACING;% [deg]  %from the center for square 4
ECC_VER5 = ECC_VER+3*SPACING;% [deg]  %from the center for square 5
ECC_VER6 = ECC_VER-3*SPACING;% [deg]  %from the center for square 6  

%-----------

% Fixation point Rect (it's red in the original stimulus)
rFixation = round(Ycenter-PixPerMin*(FIX_DIA/2)):round(Ycenter+PixPerMin*(FIX_DIA/2));
cFixation = round(Xcenter-PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX):round(Xcenter+PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX);


% Vernier gap
inner = 0.5*PixPerSec*VERNIER_GAP;       % vernier gap (divided by 2)
outer = inner+PixPerSec*VERNIER_LEN;     % inner+VERNIER LEN


% Baseline shape coordinates
xShape = [1 -1 1 1;
    1 -1 -1 -1];
yShape = [1 1 1 -1;
    -1 -1 1 -1];

% 7-star coordinates
NoPoints = 7;
Amplitude = 0.3;
phi = 2*pi;
m = (pi/2-phi)/NoPoints:pi/NoPoints:2*pi;
rm = 0.85*(1+Amplitude*sin(NoPoints*m+phi));
[x7,y7] = pol2cart(m(1:end-1),rm(1:end-1));
x7 = [x7(1:end-1); x7(2:end)];
y7 = [y7(1:end-1); y7(2:end)];


for Side = side
    
    img = zeros(height,width);
    
    % Fixation point
    img(rFixation,cFixation) = 1;
    
    %--- Draw Vernier
    vofs = (Side*Offst);                       % Vernier offset for this image
    xVern = round(vofs*0.5*[1 -1] + (Xcenter+PixPerDeg*(ECC_VER)));
    yVern = round([-outer -inner inner outer] + Ycenter);
    rVern1 = yVern(1):yVern(2);
    rVern2 = yVern(3):yVern(4);
    cVern1 = xVern(1);
    cVern2 = xVern(2);
    
    img(rVern1,cVern1) = 1;
    img(rVern2,cVern2) = 1;
    
    % Save the coordinates
    XCoords = [cVern1 cVern2; cVern1 cVern2];
    YCoords = [rVern1(1) rVern2(1); rVern1(end) rVern2(end)]; 
    
    if Cond==2 %one centered square
            
            %--- Draw a square
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; 
            yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);                          
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
    elseif Cond==3 % 7 horizontal squares
            
            %--- Draw a square
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape;
            yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            %--- Draw the other squares
            xCoords = Xcenter +PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*xShape; 
            yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*xShape; 
            yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; 
            yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; 
            yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; 
            yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; 
            yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            
    elseif Cond==4 % 7 horizontal squares and diamonds alternated
            
            
            
            %--- Draw a central square
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; 
            yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            %--- Draw the other squares
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; 
            yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; 
            yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            %Draw the stars
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; 
            yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; 
            yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; 
            yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; 
            yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; 
            yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; 
            yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; 
            yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; 
            yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            
            
    elseif Cond==5 % Pattern (3 vertical squares alternated with 3 vertical diamonds)
            
            %--- Draw the squares
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            %Draw the stars
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            
            
    elseif Cond==6 % 3 rows of 7 horizontal squares and stars alternated
            
            
            %--- Draw 1 central square
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            %Column 2
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            %Column 3
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            %Column 4
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            %Draw the 7 -stars
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            
            
    elseif Cond==7 % 3 rows of 7 horizontal squares and diamonds randomly arranged
            
            
            %--- Draw 1 central square
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img);
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            %Column 2
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            %Column 3
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            %Column 4
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            %Draw the 7 -stars
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            
            
    elseif Cond==8 % Pattern (3 vertical squares alternated with 3 vertical diamonds)
            
            %--- Draw the squares
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            %Draw the stars
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            
    elseif Cond==9 % Pattern (3 vertical squares alternated with 3 vertical diamonds)
            
            %--- Draw the squares
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            %Draw the stars
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
                     
            
    elseif Cond==10 % 7 rows of 7 horizontal squares and stars alternated
            
            %--- Draw 3 central squares
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            %Column 2
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            xCoords = Xcenter +PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            %Column 3
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            %Column 4
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            %Draw the 7 -stars
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            
    elseif Cond==11 % 7 rows of 7 horizontal squares and stars alternated
            
            %--- Draw 3 central squares
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            %Column 2
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            %Column 3
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            %Column 4
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter +PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            %Draw the 7 -stars
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER5+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER6+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER3+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER4+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
    elseif Cond==12
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            
    elseif Cond==13 % Pattern (3 vertical squares alternated with 3 vertical diamonds)
            
            %--- Draw the squares
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
            
    elseif Cond==14 % Pattern (3 vertical squares alternated with 3 vertical diamonds)
            
            %--- Draw the squares
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER1+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER3+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter +PixPerDeg*ECC_VER4+PixPerSec*INT_FLANKER*xShape; yCoords = Ycenter+PixPerDeg*ECC_VER2+PixPerSec*(FLANKER_LEN/2)*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            %Draw the stars
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER1+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER2+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER5+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            xCoords = Xcenter+PixPerDeg*ECC_VER6+PixPerSec*INT_FLANKER*x7; yCoords = Ycenter+PixPerDeg*ECC_VER+PixPerSec*(FLANKER_LEN/2)*y7;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
            
            
    end
        
    Img(:,:) = img;   
    
end