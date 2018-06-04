% Draw and save LPSY crowding stimuli. 
% The first argument indicates which stimuli to create. Use 'all' to create
% them all.
% Second argument is vernier offset. 1 is left and -1 is right.
% Type 'jpg' as the third argument to save the stimuli in .jpg, 'mat' for .mat.
%
% Usage: Drawstim('squares', 1, 'jpg');
%
% Created by Adrien Doerig 
% 7th April 2016, EPFL, LPSY Herzog Lab 

function [] = DrawStim(StimType, offset, varargin)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% directory handling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [progPath,~,~] = fileparts(which(mfilename)); % get program directory
    cd(progPath); % go there just in case we're not already

    % if there's no stimuli directory, create one
    if exist([progPath filesep '\stimuli'],'dir') == 0
        mkdir('stimuli');
    end

    stimPath = [progPath, '\stimuli'];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% draw/save everything
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch lower(StimType)
        case 'vernier'
            cd(progPath);
            [Img(:,:),~] = CtrFlnksWithVernierNoGuideMod(1,0.5, offset);
            cd(stimPath);
            figure(1)
            imshow(Img); axis image; colormap gray;
            truesize;
            if strcmpi(varargin, 'jpg')
                saveas(gcf,'vernier.jpg')
            end
            if strcmpi(varargin, 'mat')
                save('vernier.mat','Img')
            end
            if strcmpi(varargin, 'png')
                saveas(gcf,'vernier.png')
            end
        case 'malania short'
            for i = 1:4
                cd(progPath);
                [Img(:,:),~] = CtrFlnksWithVernierNoGuideMod(i,0.5, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('malaniaShort_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('malaniaShort_%d.jpg', i);
                    saveas(gcf,str);
                end
            end
        case 'malania equal'
            for i = 1:4
                cd(progPath);
                [Img(:,:),~] = CtrFlnksWithVernierNoGuideMod(i,1, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('malaniaEqual_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('malaniaEqual_%d.jpg', i);
                    saveas(gcf,str);
                end
            end      
        case 'malania long'
            for i = 1:4
                cd(progPath);
                [Img(:,:),~] = CtrFlnksWithVernierNoGuideMod(i,2, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('malaniaLong_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('malaniaLong_%d.jpg', i-1);
                    saveas(gcf,str);
                end
            end
        case 'gestalts'
            cd(progPath);
            [gestImg,~] = Gestalts(offset);
            cd(stimPath);
            for i = 1:size(gestImg,3)
                figure(1)
                imshow(gestImg(:,:,i)); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('gestalts_%d.jpg', i);
                    saveas(gcf,str);        
                end
                if strcmpi(varargin, 'mat')
                    str=sprintf('gestalts_%d.mat', i);
                    Img = gestImg(:,:,i);
                    save(str,'Img');
                end
            end  
        case 'circles'
            for i = 1:6
                cd(progPath);
                [Img(:,:),~] = MauroCirclesMod(i,offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('circles_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('circles_%d.jpg', i);
                    saveas(gcf,str);
                end
            end
        case 'squares'
            for i = 1:3
                cd(progPath);
                [Img(:,:),~] = MauroPattern2Mod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('squares_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('squares_%d.jpg', i-1);
                    saveas(gcf,str);
                end
            end
        case 'hexagons'
            for i = 1:11
                cd(progPath);
                [Img(:,:),~] = MauroHexagonsMod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    if i <10
                        str=sprintf('hexagons_0%d.mat', i);
                        save(str,'Img')   
                    else
                        str=sprintf('hexagons_%d.mat', i);
                        save(str,'Img')  
                    end
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    if i <10
                        str=sprintf('hexagons_0%d.jpg', i);
                        saveas(gcf,str);   
                    else
                        str=sprintf('hexagons_%d.jpg', i);
                        saveas(gcf,str);  
                    end
                end
            end
        case 'octagons'
            for i = 1:11
                cd(progPath);
                [Img(:,:),~] = MauroOctagonsMod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    if i <10
                        str=sprintf('octagons_0%d.mat', i);
                        save(str,'Img')   
                    else
                        str=sprintf('octagons_%d.mat', i);
                        save(str,'Img')  
                    end
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    if i <10
                        str=sprintf('octagons_0%d.jpg', i);
                        saveas(gcf,str);   
                    else
                        str=sprintf('octagons_%d.jpg', i);
                        saveas(gcf,str);  
                    end
                end
            end
        case 'stars'
            for i = 1:9
                cd(progPath);
                [Img(:,:),~] = MauroStarsMod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('stars_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('stars_%d.jpg', i);
                    saveas(gcf,str);
                end
            end
        case 'irreg 1'
            for i = 1:10
                cd(progPath);
                [Img(:,:),~] = MauroIrregular1Mod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    if i <10
                        str=sprintf('irreg1_0%d.mat', i);
                        save(str,'Img')   
                    else
                        str=sprintf('irreg1_%d.mat', i);
                        save(str,'Img')  
                    end
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    if i <10
                        str=sprintf('irreg1_0%d.jpg', i);
                        saveas(gcf,str);   
                    else
                        str=sprintf('irreg1_%d.jpg', i);
                        saveas(gcf,str);  
                    end
                end
            end
        case 'irreg 2'
            for i = 1:5
                cd(progPath);
                [Img(:,:),~] = MauroIrregular2Mod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('irreg2_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('irreg2_%d.jpg', i);
                    saveas(gcf,str);
                end
            end
        case 'patternirregular'
            for i = 1:11
                cd(progPath);
                [Img,~,~,~] = MauroPatternIrregularMod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    if i <10
                        str=sprintf('patternIrregular_0%d.mat', i);
                        save(str,'Img')   
                    else
                        str=sprintf('patternIrregular_%d.mat', i);
                        save(str,'Img')  
                    end
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    if i <10
                        str=sprintf('patternIrregular_0%d.jpg', i);
                        saveas(gcf,str);   
                    else
                        str=sprintf('patternIrregular_%d.jpg', i);
                        saveas(gcf,str);  
                    end
                end
            end
        case 'pattern 2'
            for i = 1:10
                cd(progPath);
                [Img,~,~,~] = MauroPattern2Mod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    if i <10
                        str=sprintf('pattern2_0%d.mat', i);
                        save(str,'Img')   
                    else
                        str=sprintf('pattern2_%d.mat', i);
                        save(str,'Img')  
                    end
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    if i <10
                        str=sprintf('pattern2_0%d.jpg', i);
                        saveas(gcf,str);   
                    else
                        str=sprintf('pattern2_%d.jpg', i);
                        saveas(gcf,str);  
                    end
                end
            end
        case 'pattern stars'
            for i = 1:14
                cd(progPath);
                [Img,~,~,~] = MauroPatternStarsMod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    if i <10
                        str=sprintf('patternStars_0%d.mat', i);
                        save(str,'Img')   
                    else
                        str=sprintf('patternStars_%d.mat', i);
                        save(str,'Img')  
                    end
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    if i <10
                        str=sprintf('patternStars_0%d.jpg', i);
                        saveas(gcf,str);   
                    else
                        str=sprintf('patternStars_%d.jpg', i);
                        saveas(gcf,str);  
                    end
                end
            end        

        case 'all'
            cd(progPath);
            [Img(:,:),~] = CtrFlnksWithVernierNoGuideMod(1,0.5, offset);
            cd(stimPath);
            figure(1)
            imshow(Img); axis image; colormap gray;
            truesize;
            if strcmpi(varargin, 'jpg')
                saveas(gcf,'vernier.jpg')
            end
            if strcmpi(varargin, 'mat')
                save('vernier.mat','Img')
            end
            if strcmpi(varargin, 'png')
                saveas(gcf,'vernier.png')
            end
            
            for i = 1:4
                cd(progPath);
                [Img(:,:),~] = CtrFlnksWithVernierNoGuideMod(i,0.5, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('malaniaShort_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('malaniaShort_%d.jpg', i);
                    saveas(gcf,str);
                end
            end
            
            for i = 1:4
                cd(progPath);
                [Img(:,:),~] = CtrFlnksWithVernierNoGuideMod(i,1, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('malaniaEqual_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('malaniaEqual_%d.jpg', i);
                    saveas(gcf,str);
                end
            end
            
            for i = 1:4
                cd(progPath);
                [Img(:,:),~] = CtrFlnksWithVernierNoGuideMod(i,2, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('malaniaLong_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('malaniaLong_%d.jpg', i-1);
                    saveas(gcf,str);
                end
            end
            
            cd(progPath);
            [gestImg,~] = Gestalts(offset);
            cd(stimPath);
            for i = 1:size(gestImg,3)
                figure(1)
                imshow(gestImg(:,:,i)); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('gestalts_%d.jpg', i);
                    saveas(gcf,str);        
                end
                if strcmpi(varargin, 'mat')
                    str=sprintf('gestalts_%d.mat', i);
                    Img = gestImg(:,:,i);
                    save(str,'Img');
                end
            end
            
            Img = [];
            for i = 1:6
                cd(progPath);
                [Img(:,:),~] = MauroCirclesMod(i,offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('circles_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('circles_%d.jpg', i);
                    saveas(gcf,str);
                end
            end
            
            Img = [];
            for i = 1:3
                cd(progPath);
                [Img(:,:),~] = MauroPattern2Mod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('squares_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('squares_%d.jpg', i-1);
                    saveas(gcf,str);
                end
            end
            
            Img = [];
            for i = 1:11
                cd(progPath);
                [Img(:,:),~] = MauroHexagonsMod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    if i <10
                        str=sprintf('hexagons_0%d.mat', i);
                        save(str,'Img')   
                    else
                        str=sprintf('hexagons_%d.mat', i);
                        save(str,'Img')  
                    end
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    if i <10
                        str=sprintf('hexagons_0%d.jpg', i);
                        saveas(gcf,str);   
                    else
                        str=sprintf('hexagons_%d.jpg', i);
                        saveas(gcf,str);  
                    end
                end
            end
            
            Img = [];
            for i = 1:11
                cd(progPath);
                [Img(:,:),~] = MauroOctagonsMod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    if i <10
                        str=sprintf('octagons_0%d.mat', i);
                        save(str,'Img')   
                    else
                        str=sprintf('octagons_%d.mat', i);
                        save(str,'Img')  
                    end
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    if i <10
                        str=sprintf('octagons_0%d.jpg', i);
                        saveas(gcf,str);   
                    else
                        str=sprintf('octagons_%d.jpg', i);
                        saveas(gcf,str);  
                    end
                end
            end
            
            Img = [];
            for i = 1:9
                cd(progPath);
                [Img(:,:),~] = MauroStarsMod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('stars_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('stars_%d.jpg', i);
                    saveas(gcf,str);
                end
            end
            
            Img = [];
            for i = 1:10
                cd(progPath);
                [Img(:,:),~] = MauroIrregular1Mod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    if i <10
                        str=sprintf('irreg1_0%d.mat', i);
                        save(str,'Img')   
                    else
                        str=sprintf('irreg1_%d.mat', i);
                        save(str,'Img')  
                    end
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    if i <10
                        str=sprintf('irreg1_0%d.jpg', i);
                        saveas(gcf,str);   
                    else
                        str=sprintf('irreg1_%d.jpg', i);
                        saveas(gcf,str);  
                    end
                end
            end
            
            Img = [];
            for i = 1:5
                cd(progPath);
                [Img(:,:),~] = MauroIrregular2Mod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    str=sprintf('irreg2_%d.mat', i);
                    save(str,'Img')            
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    str=sprintf('irreg2_%d.jpg', i);
                    saveas(gcf,str);
                end
            end
            
            Img = [];
            for i = 1:11
                cd(progPath);
                [Img,~,~,~] = MauroPatternIrregularMod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    if i <10
                        str=sprintf('patternIrregular_0%d.mat', i);
                        save(str,'Img')   
                    else
                        str=sprintf('patternIrregular_%d.mat', i);
                        save(str,'Img')  
                    end
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    if i <10
                        str=sprintf('patternIrregular_0%d.jpg', i);
                        saveas(gcf,str);   
                    else
                        str=sprintf('patternIrregular_%d.jpg', i);
                        saveas(gcf,str);  
                    end
                end
            end
            
            Img = [];
            for i = 1:10
                cd(progPath);
                [Img,~,~,~] = MauroPattern2Mod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    if i <10
                        str=sprintf('pattern2_0%d.mat', i);
                        save(str,'Img')   
                    else
                        str=sprintf('pattern2_%d.mat', i);
                        save(str,'Img')  
                    end
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    if i <10
                        str=sprintf('pattern2_0%d.jpg', i);
                        saveas(gcf,str);   
                    else
                        str=sprintf('pattern2_%d.jpg', i);
                        saveas(gcf,str);  
                    end
                end
            end
            
            Img = [];
            for i = 1:14
                cd(progPath);
                [Img,~,~,~] = MauroPatternStarsMod(i, offset);
                cd(stimPath);
                if strcmpi(varargin, 'mat')
                    if i <10
                        str=sprintf('patternStars_0%d.mat', i);
                        save(str,'Img')   
                    else
                        str=sprintf('patternStars_%d.mat', i);
                        save(str,'Img')  
                    end
                end
                figure(1);
                imshow(Img); axis image; colormap gray;
                truesize;
                if strcmpi(varargin, 'jpg')
                    if i <10
                        str=sprintf('patternStars_0%d.jpg', i);
                        saveas(gcf,str);   
                    else
                        str=sprintf('patternStars_%d.jpg', i);
                        saveas(gcf,str);  
                    end
                end
            end
    end

    cd(progPath);        
end
