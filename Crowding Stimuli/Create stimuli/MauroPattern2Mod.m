function [Img,Data,XCoords,YCoords] = MauroPattern2Mod(Cond,side)

% MAUROPATTERN2 Generate the stimuli from Mauro's "pattern2"
% crowding experiment.
%
% Usage: [Img,Data,XCoords,YCoords] = MauroPattern2(Condition);
%
% Conditions may be any integer from 1 to 10
%
% e.g.
%
% for Cond = 1:10
%     [Img,Data,XCoords,YCoords] = MauroPattern2(Cond);
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
% Created by Mauro Manassi (Experiment_pattern2.m)
% Modified by Aaron Clarke, 30 October 2014, EPFL, LPSY Herzog Lab to only
% produce the stimuli and not run the experiment.

data = [148.5	1973.5	874.5	1914	1153.333333	776.3333333	1240.666667	2000	1789.5	1844;
    209.5	1996.5	625.3333333	1977.5	807	1427	840.6666667	1378.5	694	1997.5;
    112.5	1998.5	645	1993.5	710	637.5	1760	1993.5	1592	283;
    111	1704	579.5	907.5	848.3333333	1306	1826.5	1628	1348	481;
    202.5	1929.5	693	1968.5	991	1707	1671	1711	1210	1993];
Data = zeros(2,length(data(1,:)));
Data(1,:) = sum(data,1)/length(data(:,1));
Data(2,:) = std(data,1)/sqrt(length(data(:,1)));

% Stimulus conversion factors
SecPerPix = 275;                    % arc seconds per pixel
SecPerDeg = 3600;                   % seconds of arc per degree of visual angle
PixPerSec = 1/SecPerPix;            % Pixels per second of arc
PixPerMin = 60/SecPerPix;           % Pixels per minute of arc
PixPerDeg = SecPerDeg/SecPerPix;    % Pixels per degree of visual angle

Offst = 1.5;                        % Vernier offset (pix)
ECC_FIX = 11;                       % [deg] %from the center
VERNIER_LEN = 2400;                 % [arcsec], one segment
VERNIER_GAP = 240;                  % [arcsec]
ECC_VER = 0;
FLANKER_LEN = 7200;                 % [arcsec]
FIX_DIA = 8;                        % [arcmin]
sdLine = 0.5;                       % standard error of the gaussian,
                                    % corresponds to line width
INT_FLANKER = 3600;                 % [arcsec]


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

SPACING = 2.2;                  % [deg]  %from the center for the vernier
ECC_VER1 = ECC_VER+SPACING;     % [deg]  %from the center for square 1
ECC_VER2 = ECC_VER-SPACING;     % [deg]  %from the center for square 2
ECC_VER3 = ECC_VER+2*SPACING;   % [deg]  %from the center for square 3
ECC_VER4 = ECC_VER-2*SPACING;   % [deg]  %from the center for square 4
ECC_VER5 = ECC_VER+3*SPACING;   % [deg]  %from the center for square 5
ECC_VER6 = ECC_VER-3*SPACING;   % [deg]  %from the center for square 6

%-----------

% Fixation point Rect (it's red in the original stimulus)
rFixation = round(Ycenter-PixPerMin*(FIX_DIA/2)):round(Ycenter+PixPerMin*(FIX_DIA/2));
cFixation = round(Xcenter-PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX):round(Xcenter+PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX);


% Vernier gap
inner = 0.5*PixPerSec*VERNIER_GAP;       % vernier gap (divided by 2)
outer = inner+PixPerSec*VERNIER_LEN;     % inner+VERNIER LEN


% Baseline shape coordinates

% square coordinates
xSquare = [1 -1 1 1;
    1 -1 -1 -1];
ySquare = [1 1 1 -1;
    -1 -1 1 -1];

% 7-star coordinates
NoPoints = 7;
Amplitude = 0.3;
phi = 2*pi;
m = (pi/2-phi)/NoPoints:pi/NoPoints:2*pi;
rm = 0.85*(1+Amplitude*sin(NoPoints*m+phi));
[x7,y7] = pol2cart(m(1:end-1),rm(1:end-1));
xStar = [x7(1:end-1); x7(2:end)];
yStar = [y7(1:end-1); y7(2:end)];

% Draw circles
EccCirc = [-2.97 -1.98 -0.99 0 0.99 1.98 2.97]*SPACING/2.8*PixPerSec*10200;
nEccCirc = length(EccCirc);

% Set up a Cartesian coordinate system over the image
[x,y] = meshgrid(1:width,1:height);

% Circle parameters
Rad = PixPerSec*INT_FLANKER;

% Make the circles
nEllipsePoints = 50;
d = zeros(height,width,3,7);
cX = zeros(2,2*nEllipsePoints-1,3,7);
cY = cX;
for yEccCirc = 3:5
    for xEccCirc = 1:nEccCirc    
        rad = sqrt((x-Xcenter+EccCirc(xEccCirc)).^2 + ((y-Ycenter)+EccCirc(yEccCirc)).^2);
        d(:,:,yEccCirc-2,xEccCirc) = exp(-(rad-Rad).^2/(2*sdLine^2)); 
        [cx,cy] = ellipse(nEllipsePoints,Rad,Rad,Xcenter+EccCirc(xEccCirc),...
            Ycenter+EccCirc(yEccCirc),0); 
        cX(:,:,yEccCirc-2,xEccCirc) = [cx(1:end-1); cx(2:end)];
        cY(:,:,yEccCirc-2,xEccCirc) = [cy(1:end-1); cy(2:end)];
    end
end



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
    
    %--- Choose one of the conditions
    if Cond==2 %one centered square
        
        %--- Draw a square
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare;
        yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        
        
        
    elseif Cond==3 % 7 horizontal squares
        
        %--- Draw a square
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare;
        yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        %--- Draw the other squares
        xCoords = Xcenter+PixPerDeg*(ECC_VER1)+PixPerSec*(INT_FLANKER)*xSquare;
        yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER2)+PixPerSec*(INT_FLANKER)*xSquare;
        yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare;
        yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square3
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare;
        yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square4
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER5)+PixPerSec*(INT_FLANKER)*xSquare;
        yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square5
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER6)+PixPerSec*(INT_FLANKER)*xSquare;
        yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square6
        
        
    elseif Cond==4 % 7 horizontal squares and diamonds alternated
        
        % Concatenate the circles together        
        img = cat(3,reshape(d(:,:,2,[1 3 5 7]),[height width 4]),img);
        XCoords = cat(2,XCoords,cX(:,:,2,1),cX(:,:,2,3),cX(:,:,2,5),cX(:,:,2,7));
        YCoords = cat(2,YCoords,cY(:,:,2,1),cY(:,:,2,3),cY(:,:,2,5),cY(:,:,2,7));
        
        % Take the max
        img = max(img,[],3);
        
        %--- Draw a central square
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); 
        XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        %--- Draw the other squares
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); 
        XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); 
        XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        
        
    elseif Cond==5 % Pattern (3 vertical squares alternated with 3 vertical diamonds)
        
        % Concatenate the circles together
        img = cat(3,d(:,:,1,1),d(:,:,1,3),d(:,:,1,5),d(:,:,1,7),...
            d(:,:,2,1),d(:,:,2,3),d(:,:,2,5),d(:,:,2,7),...
            d(:,:,3,1),d(:,:,3,3),d(:,:,3,5),d(:,:,3,7),img);
        XCoords = cat(2,XCoords,cX(:,:,1,1),cX(:,:,1,3),cX(:,:,1,5),cX(:,:,1,7),...
            cX(:,:,2,1),cX(:,:,2,3),cX(:,:,2,5),cX(:,:,2,7),...
            cX(:,:,3,1),cX(:,:,3,3),cX(:,:,3,5),cX(:,:,3,7));
        YCoords = cat(2,YCoords,cY(:,:,1,1),cY(:,:,1,3),cY(:,:,1,5),cY(:,:,1,7),...
            cY(:,:,2,1),cY(:,:,2,3),cY(:,:,2,5),cY(:,:,2,7),...
            cY(:,:,3,1),cY(:,:,3,3),cY(:,:,3,5),cY(:,:,3,7));
        
        % Take the max
        img = max(img,[],3);
        
        %--- Draw 3 vertical squares, central column
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); 
        XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); 
        XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); 
        XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        %--- Draw 3 vertical squares, 3rd column
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        %--- Draw 3 vertical squares, 3rd column
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        
    elseif Cond==6 % 3 rows of 7 horizontal squares and diamonds alternated
        
        
        % Concatenate the circles together
        img = cat(3,d(:,:,2,5),d(:,:,2,3),d(:,:,2,7),d(:,:,2,1),d(:,:,3,4),...
            d(:,:,3,6),d(:,:,3,2),d(:,:,1,4),d(:,:,1,6),d(:,:,1,2),img);
        XCoords = cat(2,XCoords,cX(:,:,2,5),cX(:,:,2,3),cX(:,:,2,7),cX(:,:,2,1),cX(:,:,3,4),...
            cX(:,:,3,6),cX(:,:,3,2),cX(:,:,1,4),cX(:,:,1,6),cX(:,:,1,2));
        YCoords = cat(2,YCoords,cY(:,:,2,5),cY(:,:,2,3),cY(:,:,2,7),cY(:,:,2,1),cY(:,:,3,4),...
            cY(:,:,3,6),cY(:,:,3,2),cY(:,:,1,4),cY(:,:,1,6),cY(:,:,1,2));
        
        % Take the max
        img = max(img,[],3);
        
        %--- Draw 1 central square
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        %Column 2
        xCoords = Xcenter+PixPerDeg*(ECC_VER1)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER1)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER2)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER2)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        %Column 3
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        %Column 4
        xCoords = Xcenter+PixPerDeg*(ECC_VER5)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER5)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER6)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER6)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        
        
    elseif Cond==7 % 3 rows of 7 horizontal squares and diamonds randomly arranged
        
        % Concatenate the circles together
        img = cat(3,d(:,:,2,5),d(:,:,2,3),d(:,:,2,7),d(:,:,2,1),d(:,:,1,4),...
            d(:,:,1,2),d(:,:,3,3),d(:,:,3,4),d(:,:,3,6),d(:,:,1,5),img);
        XCoords = cat(2,XCoords,cX(:,:,2,5),cX(:,:,2,3),cX(:,:,2,7),cX(:,:,2,1),cX(:,:,1,4),...
            cX(:,:,1,2),cX(:,:,3,3),cX(:,:,3,4),cX(:,:,3,6),cX(:,:,1,5));
        YCoords = cat(2,YCoords,cY(:,:,2,5),cY(:,:,2,3),cY(:,:,2,7),cY(:,:,2,1),cY(:,:,1,4),...
            cY(:,:,1,2),cY(:,:,3,3),cY(:,:,3,4),cY(:,:,3,6),cY(:,:,1,5));
        
        % Take the max
        img = max(img,[],3);
        
        %--- Draw 1 central square
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        %Column 2
        xCoords = Xcenter+PixPerDeg*(ECC_VER1)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER2)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        %Column 3
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        %Column 4
        xCoords = Xcenter+PixPerDeg*(ECC_VER5)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER5)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER6)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER6)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        
        
    elseif Cond==8 % Pattern (3 vertical squares alternated with 3 vertical diamonds)
        
        % Concatenate the circles together
        img = cat(3,d(:,:,2,5),d(:,:,2,3),d(:,:,2,7),d(:,:,2,1),d(:,:,3,5),...
            d(:,:,3,3),d(:,:,3,7),d(:,:,3,1),d(:,:,1,5),d(:,:,1,3),d(:,:,1,7),d(:,:,1,1),img);
        XCoords = cat(2,XCoords,cX(:,:,2,5),cX(:,:,2,3),cX(:,:,2,7),cX(:,:,2,1),cX(:,:,3,5),...
            cX(:,:,3,3),cX(:,:,3,7),cX(:,:,3,1),cX(:,:,1,5),cX(:,:,1,3),cX(:,:,1,7),cX(:,:,1,1));
        YCoords = cat(2,YCoords,cY(:,:,2,5),cY(:,:,2,3),cY(:,:,2,7),cY(:,:,2,1),cY(:,:,3,5),...
            cY(:,:,3,3),cY(:,:,3,7),cY(:,:,3,1),cY(:,:,1,5),cY(:,:,1,3),cY(:,:,1,7),cY(:,:,1,1));
        
        % Take the max
        img = max(img,[],3);
        
        %--- Draw 3 vertical squares, central column
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        %--- Draw 3 vertical squares, 3rd column
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        %--- Draw 3 vertical squares, 3rd column
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        
        
        
    elseif Cond==9 % Pattern (3 vertical squares alternated with 3 vertical circles)+ 1 star in position (1,1)
        
        
        % Concatenate the circles together
        img = cat(3,d(:,:,2,5),d(:,:,2,3),d(:,:,2,7),d(:,:,2,1),d(:,:,3,4),...
            d(:,:,3,6),d(:,:,3,2),d(:,:,1,4),d(:,:,1,6),d(:,:,1,2),img);
        XCoords = cat(2,XCoords,cX(:,:,2,5),cX(:,:,2,3),cX(:,:,2,7),cX(:,:,2,1),cX(:,:,3,4),...
            cX(:,:,3,6),cX(:,:,3,2),cX(:,:,1,4),cX(:,:,1,6),cX(:,:,1,2));
        YCoords = cat(2,YCoords,cY(:,:,2,5),cY(:,:,2,3),cY(:,:,2,7),cY(:,:,2,1),cY(:,:,3,4),...
            cY(:,:,3,6),cY(:,:,3,2),cY(:,:,1,4),cY(:,:,1,6),cY(:,:,1,2));
        
        % Take the max
        img = max(img,[],3);
        
        
        %--- Draw 1 central square
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER1)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER2)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER2)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        %Column 3
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        %Column 4
        xCoords = Xcenter+PixPerDeg*(ECC_VER5)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER5)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER6)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER6)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        %Draw a 7-star
        xCoords = Xcenter+PixPerDeg*(ECC_VER1)+PixPerSec*(INT_FLANKER)*xStar;
        yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*yStar;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER1)+PixPerSec*(INT_FLANKER)*xStar;
        yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*yStar;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        
    elseif Cond==10 %7 alternated squares and circles + 3 vertical squares
        
        % Concatenate the circles together
        img = cat(3,d(:,:,2,5),d(:,:,2,3),d(:,:,2,7),d(:,:,2,1),img);
        XCoords = cat(2,XCoords,cX(:,:,2,5),cX(:,:,2,3),cX(:,:,2,7),cX(:,:,2,1));
        YCoords = cat(2,YCoords,cY(:,:,2,5),cY(:,:,2,3),cY(:,:,2,7),cY(:,:,2,1));
        
        % Take the max
        img = max(img,[],3);
        
        
        %--- Draw 3 vertical squares, central column
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER1)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square1
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerDeg*(ECC_VER2)+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords); %square2
        
        %Draw the other squares
        xCoords = Xcenter+PixPerDeg*(ECC_VER3)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter+PixPerDeg*(ECC_VER4)+PixPerSec*(INT_FLANKER)*xSquare; yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*ySquare;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine); img = max(ShapeImg,img); XCoords = cat(2,XCoords,xCoords); YCoords = cat(2,YCoords,yCoords);
        
        
    end
    Img(:,:) = img;    
end