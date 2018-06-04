% a function to run the model with a specific stimulus (given as the argument)
% creates a results folder to save the output figure.
%
% Created by Adrien Doerig, LPSY, september 2016

function [CrossCorr] = aaronModel(stim)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % import needed variables from the base workspace
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fittingStim = evalin('base','fittingStims{run}'); % !!!NEED DIFFERENT NAMES IF RUNNING RUN ALL OR FITGLOBALPSYCHFUN get fittingStim from the base workspace into the function
    FitParams = evalin('base', 'FitParams');
    fitting = evalin('base', 'fitting');

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % directory management
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [progPath,~,~] = fileparts(which(mfilename)); % get program directory
    cd(progPath); % go there just in case we're not already

    % if there's no results directory, it's a new stimulus
    if exist([progPath filesep '\results\psychometric function fitted on ' fittingStim],'dir') == 0
        resDir = 0;

    % if a directory already exists, then it's a returning stimulus
    elseif exist([progPath filesep '\results\psychometric function fitted on ' fittingStim],'dir') == 7
        resDir = 1;
    end

    if resDir == 0
        mkdir('results');
        cd([progPath, '\results'])
        mkdir(['psychometric function fitted on ' fittingStim])
    end

    % define paths
    resPath = [progPath, '\results\psychometric function fitted on ' fittingStim]; % path to data folder
    stimPath = [progPath, '\stimuli'];

    % create filename that iterates to avoid overwriting
    cd(resPath); % go to data path
    d=dir('*.mat'); % isolate all .mat files (i.e. data files)
    numfiles=length(d);
    gotName=0;
    N=0;
    while(gotName==0)
        N=N+1;
        % set up desired filename
        fileName = [stim,'_',num2str(N)];
        foundName=0;
        for kk=1:numfiles % loop through every mat file and look for the currently desired filename
            tmp=d(kk).name;
            foundName=strcmpi(filename,tmp);
            if(foundName==1)
                break % if we find it, go back and iterate
            end
        end
        if foundName==0
            gotName=1; % if we don't find it, use this filename
        end
    end
    cd(progPath); % go back to root

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % setup stimuli, data, etc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd(stimPath);
    files = dir([stimPath, '\', stim, '*.mat']);
    names = {files.name};
    Img = [];
    for i = 1:length(names)
        img = load(names{i});
        Img(:,:,i) = img.Img;
    end
    nCond = size(Img,3); 
    
    cd(progPath);
    
    if ~fitting
        makeData
        Data = dataStruct.(matlab.lang.makeValidName(stim));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run the model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    endStoppedModel

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~fitting % no need to plot, etc. if we are fitting the psych function
        cd(resPath);
%         saveas(gcf,[stim, '.fig'])
        saveas(gcf,[stim, '.jpg'])
        cd(progPath);
    end

end