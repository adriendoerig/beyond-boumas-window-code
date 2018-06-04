function [Img,Data,XCoords,YCoords] = MauroPatternIrregular(Cond,side)

% MAUROPATTERNIRREGULAR Generate the stimuli from Mauro's "pattern_irregular"
% crowding experiment.
%
% Usage: [Img,Data,XCoords,YCoords] = MauroPatternIrregular(Condition);
%
% Conditions may be any integer from 1 to 11
%
% e.g.
%
% for Cond = 1:11
%     [Img,Data,XCoords,YCoords] = MauroPatternIrregular(Cond);
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
% Created by Mauro Manassi (Experiment_pattern_irregular.m)
% Modified by Aaron Clarke, 24 October 2014, EPFL, LPSY Herzog Lab to only
% produce the stimuli and not run the experiment.

% Human Data
data = [98	1870	559	437.5	1842	1870	1969.5	1995	1670.666667	1750	1860;
    262	1488.666667	1130.333333	1077.5	1299.5	503.5	539	1365	1374	1344.666667	1622.666667;
    121	1760	1348	802.6666667	267.5	732.3333333	935	384.6666667	688.6666667	1029	756;
    59	1953.5	123.5	176.5	453	88.5	75	59	252	400.5	162;
    134	1692	1151.75	1299	1039	1754.5	1447	1360.5	1171.333333	1029	1064.333333;
    135.5	1769.5	1194	1309	1109.5	1592.5	1317	815.3333333	1396.666667	1767.5	1756;
    82	1857.333333	1057	1303	1551.333333	1397	1290.5	1246.333333	1956.5	1931	1439.5];
Data = zeros(2,length(data(1,:)));
Data(1,:) = sum(data,1)/length(data(:,1));
Data(2,:) = std(data,1)/sqrt(length(data(:,1)));

clarkeFactor = 3.56; % match expected size in endstoppedmodel

% Stimulus conversion factors
SecPerPix = 275;                    % arc seconds per pixel
SecPerDeg = 3600;                   % seconds of arc per degree of visual angle
PixPerSec = 1/SecPerPix;            % Pixels per second of arc
PixPerMin = 60/SecPerPix;           % Pixels per minute of arc
PixPerDeg = SecPerDeg/SecPerPix;    % Pixels per degree of visual angle

Offst = 1.5*2;                        % Vernier offset (pix)
ECC_FIX = 9*clarkeFactor;                        % [deg] %from the center
VERNIER_LEN = 2400*clarkeFactor;                 % [arcsec], one segment
VERNIER_GAP = 240*clarkeFactor;                  % [arcsec]
POINTER_GAP = 18000*clarkeFactor;                % [arcsec] for the pointers
ECC_VER = 0*clarkeFactor;
FLANKER_LEN = 7200*clarkeFactor;                 % [arcsec]
FIX_DIA = 8*clarkeFactor;                        % [arcmin]
sdLine = 0.5;                       % standard error of the gaussian,
                                    % corresponds to line width
IRR_R = 1.4*(FLANKER_LEN/2);

%-----------

% Image dimensions
height = 500;
width = 1000;
ImSz = [height width];
Img = zeros(height,width,1);

Xcenter = 0.5*width;   % center X coordinates
Ycenter = 0.5*height;  % center Y coordinates

SPACING = 2.7*clarkeFactor;                  % [deg]  %from the center for the vernier
Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing

%-----------

% Fixation point Rect (it's red in the original stimulus)
rFixation = round(Ycenter-PixPerMin*(FIX_DIA/2)):round(Ycenter+PixPerMin*(FIX_DIA/2));
cFixation = round(Xcenter-PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX):round(Xcenter+PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX);


% Vernier gap
inner = 0.5*PixPerSec*VERNIER_GAP;       % vernier gap (divided by 2)
outer = inner+PixPerSec*VERNIER_LEN;     % inner+VERNIER LEN


% Baseline shape coordinates
xShape = [-0.8 -0.4 -0.7 -0.2 0.0 0.3 0.5 0.8 0.5 0.6 0.3 0.0 -0.3 -0.6;
    -0.4 -0.7 -0.2 0.0 0.3 0.5 0.8 0.5 0.6 0.3 0.0 -0.3 -0.6 -0.8];
yShape = [-0.8 0.5 0.3 0.9 0.6 0.7 0.1 0.0 -0.5 -0.7 -0.9 -0.6 -0.9 -0.7;
    0.5 0.3 0.9 0.6 0.7 0.1 0.0 -0.5 -0.7 -0.9 -0.6 -0.9 -0.7 -0.8];
XShape = [xShape(:)'; yShape(:)'];
RotCoords = [cos(pi)  -sin(pi); sin(pi) cos(pi)]*XShape;
xShape180 = reshape(RotCoords(1,:),size(xShape));
yShape180 = reshape(RotCoords(2,:),size(yShape));


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
    
    if Cond==2 %Irregular form (100%)
        
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
    elseif Cond==3 %7 irregular shapes
        
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        % Make the flanking shapes
        for iEcc = 1:length(Ecc_Ver)
            
            xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*IRR_R*xShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
    elseif Cond==4 %7 shapes touching; 3 lines
        
        
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        % Make the flanking shapes
        for xEcc = 1:length(Ecc_Ver)
            
            xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
        end
        
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*(IRR_R)*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        
        for xEcc = 1:length(Ecc_Ver)
            
            xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*(IRR_R)*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        for xEcc = 1:length(Ecc_Ver)
            
            xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
        
    elseif Cond==5 %Flankers alternated (180 or 0 degrees)
        
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        % Make the flanking shapes
        for xEcc = 1:length(Ecc_Ver)
            
            if iseven(Ecc_Ver(xEcc)/SPACING)
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
                yCoords = Ycenter+PixPerSec*IRR_R*yShape;
            else
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape180;
                yCoords = Ycenter+PixPerSec*IRR_R*yShape180;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        
        end
        
        
    elseif Cond==6 % Pattern (3 vertical irregular2 alternated with 3 vertical irregular2_180)
        
        % 3 central vertical irregular2
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        % Make the flanking shapes
        for xEcc = 1:length(Ecc_Ver)
            
            if iseven(Ecc_Ver(xEcc)/SPACING)
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
                yCoords = Ycenter+PixPerSec*IRR_R*yShape;
            else
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape180;
                yCoords = Ycenter+PixPerSec*IRR_R*yShape180;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        % Make the flanking shapes
        for xEcc = 1:length(Ecc_Ver)
            
            if iseven(Ecc_Ver(xEcc)/SPACING)
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
                yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape;
            else
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape180;
                yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape180;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        % Make the flanking shapes
        for xEcc = 1:length(Ecc_Ver)
            
            if iseven(Ecc_Ver(xEcc)/SPACING)
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
                yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape;
            else
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape180;
                yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape180;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
        
    elseif Cond==7 %Pattern alternated; 3 lines of irregular/irregular180 alternated, shifted on every line
        
        IRR_R = 1.4*(FLANKER_LEN/2);
        
        %3 central vertical irregular2
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        % Make the flanking shapes
        for xEcc = 1:length(Ecc_Ver)
            
            if iseven(Ecc_Ver(xEcc)/SPACING)
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
                yCoords = Ycenter+PixPerSec*IRR_R*yShape;
            else
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape180;
                yCoords = Ycenter+PixPerSec*IRR_R*yShape180;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
        xCoords = Xcenter + PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape180;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape180;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        % Make the flanking shapes
        for xEcc = 1:length(Ecc_Ver)
            
            if iseven(Ecc_Ver(xEcc)/SPACING)
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape180;
                yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape180;
            else
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
                yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        
        end
        
        
        xCoords = Xcenter + PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape180;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape180;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        % Make the flanking shapes
        for xEcc = 1:length(Ecc_Ver)
            
            if iseven(Ecc_Ver(xEcc)/SPACING)
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape180;
                yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape180;
            else
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
                yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
        
    elseif Cond==8 %Pattern alternated; 3 lines of irregular/irregular180 alternated, shifted on every line + swith at position (1;1) and (1;-1)
        
        IRR_R = 1.4*(FLANKER_LEN/2);
        
        %3 central vertical irregular2
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        % Make the flanking shapes
        for xEcc = 1:length(Ecc_Ver)
            
            if iseven(Ecc_Ver(xEcc)/SPACING)
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
                yCoords = Ycenter+PixPerSec*IRR_R*yShape;
            else
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape180;
                yCoords = Ycenter+PixPerSec*IRR_R*yShape180;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        
        end
        
        xCoords = Xcenter + PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape180;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape180;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape180;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape180;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        %6vertical irregular2_180 column 2, left and right
        
        xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*xShape180;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape180;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*xShape180;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape180;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        %6vertical irregular2 column 3, left and riaght
        
        xCoords = Xcenter + 0.85*PixPerDeg*(Ecc_Ver(5))+PixPerSec*IRR_R*xShape180;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape180;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + 0.85*PixPerDeg*(Ecc_Ver(5))+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + 0.85*PixPerDeg*(Ecc_Ver(2))+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + 0.85*PixPerDeg*(Ecc_Ver(2))+PixPerSec*IRR_R*xShape180;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape180;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        %6vertical irregular2_180 column 4, left and right
        
        xCoords = Xcenter + 0.85*PixPerDeg*(Ecc_Ver(6))+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + 0.85*PixPerDeg*(Ecc_Ver(6))+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + 0.85*PixPerDeg*(Ecc_Ver(1))+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + 0.85*PixPerDeg*(Ecc_Ver(1))+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        
        
    elseif Cond==9 % Pattern (3 vertical irregular2 alternated with 3 vertical irregular2_180); without the central column
        
        IRR_R = 1.4*(FLANKER_LEN/2);
        
        %3 central vertical irregular2
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        % Make the flanking shapes
        yOffst = [0 PixPerDeg*Ecc_Ver(3) PixPerDeg*Ecc_Ver(4)];
        for xEcc = 1:length(Ecc_Ver)
            
            for iOffst = 1:length(yOffst)
                
                if iseven(Ecc_Ver(xEcc)/SPACING)
                    xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
                    yCoords = Ycenter+yOffst(iOffst)+PixPerSec*IRR_R*yShape;
                else
                    xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape180;
                    yCoords = Ycenter+yOffst(iOffst)+PixPerSec*IRR_R*yShape180;
                end
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end
            
        end
        
    elseif Cond==10 % 3 vertical irregular2 alternated with 3 vertical irregular2_180 in one line + one central colum of irregular 2
        
        IRR_R = 1.4*(FLANKER_LEN/2);
        
        %3 central vertical irregular2
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        for xEcc = 1:length(Ecc_Ver)
            
            if iseven(Ecc_Ver(xEcc)/SPACING)
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
                yCoords = Ycenter+PixPerSec*IRR_R*yShape;
            else
                xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape180;
                yCoords = Ycenter+PixPerSec*IRR_R*yShape180;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
        xCoords = Xcenter + PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
    elseif Cond==11 % E with more lines (5)
        
        IRR_R = 1.4*(FLANKER_LEN/2);
        
        %5 central vertical irregular2
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        % Make the flanking shapes
        yOffst = PixPerDeg*[0 Ecc_Ver(2:5)];
        for xEcc = 1:length(Ecc_Ver)
            
            for iOffst = 1:length(yOffst)
                
                if iseven(Ecc_Ver(xEcc)/SPACING)
                    xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape;
                    yCoords = Ycenter+yOffst(iOffst)+PixPerSec*IRR_R*yShape;
                else
                    xCoords = Xcenter + 0.85*PixPerDeg*Ecc_Ver(xEcc) + PixPerSec*IRR_R*xShape180;
                    yCoords = Ycenter+yOffst(iOffst)+PixPerSec*IRR_R*yShape180;
                end
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end
            
        end
        
        
        xCoords = Xcenter + PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(2)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(3)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(4)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        xCoords = Xcenter + PixPerDeg*ECC_VER+PixPerSec*IRR_R*xShape;
        yCoords = Ycenter+PixPerDeg*Ecc_Ver(5)+PixPerSec*IRR_R*yShape;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
    end
    
    Img(:,:) = img;
    
end