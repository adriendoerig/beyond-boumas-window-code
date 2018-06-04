% LnDetSpatAssocFldGeneral.m
% 
% End-stopped line detectors with spatial association fields.
% General purpose code for simulations.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

video = 1;
save_output_image = 1;

% adjust end-stopped RF size to stimulus size
if contains(stim,'malania') ||  contains(stim,'circles') ||  contains(stim,'gestalts')
    clarkeFactor = 3.56; 
else
    clarkeFactor = 1;
end

% Hermens model parameters ---
FinalTime = 40;             % Amount of time to simulate (ms)
dt = 0.75;                  % Time step size
nT = round(FinalTime/dt);   % Number of time steps
se = 7.5;                   % Excitation field standard deviation
si = 12.5;                  % Inhibition field standard deviation 
WeR = 31;
WiR = 46;
weBtw = 1.1;                % Keep this one fixed (it doesn't change that much anyway)
we = 0.01;
wi = 1;                     % Higher makes the curves go downward
wBtw = 0.975;               % Weight on neighbouring line lenghts - increasing it puts conditions in right order, but makes equal get better with flanker number

% other params
te = 16;            % Excitatory time constant
ti = 4;             % Inhibitory time constant

% Receptive fields/filters
conste = we/sqrt(2*pi*se*se);
consti = wi/sqrt(2*pi*si*si);


% excitatory
X = -WeR+1:WeR-1;
[xx,yy] = meshgrid(X,X);
[theta,r] = cart2pol(xx,yy);
ThMsk = ones(size(theta));
ThThresh = pi/6;
ThMsk(abs(theta)<ThThresh) = 0;
ThMsk(abs(theta)>pi-ThThresh) = 0;
ThMsk(WeR,WeR) = 1;
We = exp(-X.^2/(2*se^2));
We = conste * (We'*We).*ThMsk;

% inhibitory
X = -WiR+1:WiR-1;
[xx,yy] = meshgrid(X,X);
[theta,r] = cart2pol(xx,yy);
ThMsk = zeros(size(theta));
ThMsk(abs(theta)<ThThresh) = 1;
ThMsk(abs(theta)>pi-ThThresh) = 1;
Wi = exp(-X.^2/(2*si^2));
Wi = consti * (Wi'*Wi).*ThMsk;

% Hubel and Weisel style end-stopped V1 filter
SecPerDeg = 3600;
SecPerPix = 275;                             % arc seconds per pixel
PixPerSec = 1/SecPerPix;                    % Pixels per arc second
uFlnkHght = round([2400 4800 9600]*PixPerSec*clarkeFactor);        % Flanker Height, in arcsecs in brackets, transformed in pixels thereafter
nFlankHeight = length(uFlnkHght);
V = cell(nFlankHeight,1);

% Create the Vernier template
template = Img(:,:,1);
[r,c] = size(template);
ImSz = [r c];

% Make the line detector filters
for iFlnk = 1:nFlankHeight   
            
    V{iFlnk} = -ones(uFlnkHght(iFlnk)*2,1);
    if iseven(uFlnkHght(iFlnk))
        Hlf = uFlnkHght(iFlnk)/2;
        V{iFlnk}(uFlnkHght(iFlnk)-Hlf+1:uFlnkHght(iFlnk)+Hlf,1) = 1;
    else
        Hlf = (uFlnkHght(iFlnk)-1)/2;
        V{iFlnk}(uFlnkHght(iFlnk)-Hlf:uFlnkHght(iFlnk)+Hlf,1) = 1;
    end    
    
    % Normalize by number of elements
    V{iFlnk} = V{iFlnk}/(uFlnkHght(iFlnk)*2);
    
end

% Between filter-size inhibitory field
WiBtw = zeros(nFlankHeight,nFlankHeight);
sBtw = 30;
for iFilt = 1:nFlankHeight
    WiBtw(iFilt,:) = wBtw*exp(-(uFlnkHght-uFlnkHght(iFilt)).^2/sBtw^2);    
end
WiBtw = reshape(WiBtw',[1 1 nFlankHeight nFlankHeight]);
WiBtw = repmat(WiBtw,[ImSz 1 1]);

% Initialize the excitatory activations matrix
AllAe = zeros([ImSz nFlankHeight nCond]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% network dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iCond = 1:nCond

    % Generate the stimulus    
    img = Img(:,:,iCond)*2-1;
    
    % Initialize the input matrix
    inp = zeros([ImSz nFlankHeight]);
    
    % Convolve the image with the filters
    for iFilt = 1:nFlankHeight   
        inp(:,:,iFilt) = ZeroPad(conv2(img,V{iFilt},'same'),ImSz);        
    end

    % Initialize the excitatory and inhibitory activations matricies
    Ae = zeros([ImSz nFlankHeight]);
    Ai = zeros([ImSz nFlankHeight]);
    
    % Initialize the weight matricies
    wee = Ae;
    wii = Ae;        
    
    if video || save_output_image
        figure(2000)
        subplot(4,3,11)
        imagesc(img)
        title('Input')
    end
    
    % Implement dynamics
    for t = 1:nT

        % Perform the within-filter size convolutions
        for iFilt = 1:nFlankHeight
            wee(:,:,iFilt) = 0.5*conv2(Ae(:,:,iFilt),We,'same');
            wii(:,:,iFilt) = -0.5*conv2(Ai(:,:,iFilt),Wi,'same');
        end
        
        % Perform the between-filter size convolutions
        for iFilt = 1:nFlankHeight
            wee(:,:,iFilt) = sum(weBtw*wee.*WiBtw(:,:,:,iFilt),3);
            wii(:,:,iFilt) = sum(wii.*WiBtw(:,:,:,iFilt),3);
        end
            
        
        Ae = Ae + dt*(-Ae + he(wee+wii+inp))/te;
        Ai = Ai + dt*(-Ai + hi(wee+wii+inp))/ti;
        
        if video || save_output_image
            figure(2000)
            for i = 1:3
                subplot(4,3,i)
                imagesc(Ae(:,:,i));
                title(['Excitatory ', num2str(i)])
                subplot(4,3,i+3)
                imagesc(Ai(:,:,i));
                title(['Inhibitory ', num2str(i)])
                drawnow
                if iCond>1
                    subplot(4,3,i+6)
                    imagesc(AllAe(:,:,i,1));
                    xcorr = sum(sum(Ae(:,:,i).*AllAe(:,:,i,1)))/max(1,sum(sum(AllAe(:,:,i,1))));
                    title(['Vernier alone. CrossCorr = ', num2str(xcorr)])
                    axis image
                    colormap hot
                    colorbar
                end
            end
        end
    end
    
    if save_output_image
        cd(resPath);
%         saveas(gcf,[stim, '.fig'])
        saveas(gcf,[stim, '_', num2str(iCond), '_output_image.jpg'])
        cd(progPath);
    end

    AllAe(:,:,:,iCond) = Ae;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% relationship to crowding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the cross correlations of the vernier with the stimuli
CrossCorr = zeros(1,nCond);
for iCond = 1:nCond
    CrossCorr(iCond) = sum(sum(sum(AllAe(:,:,:,iCond).*AllAe(:,:,:,1))));
%     for iFlk = 1:3
%         CrossCorr(iCond) = CrossCorr(iCond) + sum(sum(AllAe(:,:,iFlk,iCond).*AllAe(:,:,iFlk,1)))/sum(sum(AllAe(:,:,iFlk,1)));
%     end
end

% Model threshold elevations 
CrossCorr = CrossCorr(2:end)/CrossCorr(1);

if ~fitting % no need to plot, etc. if we are fitting the psych function
    P = SigmoidML(FitParams,CrossCorr);
    Residuals = sqrt((P-Data(1,:)).^2);

    % plot
    figure(1);
    clf
    histogramStimsOutput = SigmoidML(FitParams,1);
    histogramStimsOutput = [histogramStimsOutput, P];
    histogramThreshElev = Data(1,:);

    bar([[110 histogramThreshElev].', histogramStimsOutput.']); % the first number is the human data vernier threshold
    legend('humans', 'Wilson & Cowan network')
    title([stim, ' with ', fittingStim, ' psychometric function'])
    
%     subplot(3,1,1);
%     BarH = bar(1:size(Data,2),Data(1,:)); set(BarH,'FaceColor','k');
%     hold on; errorbar(1:size(Data,2),Data(1,:),Data(2,:),'k.'); hold off;
%     title([stim,': Sum of Residuals = ', num2str(sum(Residuals))]);
% 
%     subplot(3,1,2);
%     BarH = bar(1:size(Data,2),P); set(BarH,'FaceColor','k');
% 
%     for ind = 1:size(Img,3)-1
%         subplot(3,size(Img,3)-1,ind+2*(size(Img,3)-1));
%         imagesc(Img(:,:,ind+1));
%         axis image;
%         axis off;
%     end
%     colormap gray;
end

