% use this to run all stimuli through the model
%
% A.Doerig, LPSY, September 2016

clear 
close all

% stims = {'circles', 'gestalts', 'hexagons', 'irreg1', 'irreg2', ...
%     'malaniaEqual', 'malaniaLong', 'malaniaShort', 'octagons' ...
%     'pattern2', 'patternIrregular', 'patternStars', 'squares', 'stars'};
stims = {'squares', 'stars'};
fittingStims = {'squares'};
% fittingStims = {'circles', 'gestalts', 'hexagons', 'irreg1', 'irreg2', ...
%     'malaniaEqual', 'malaniaLong', 'malaniaShort', 'malania', 'octagons' ...
%     'pattern2', 'patternIrregular', 'patternStars', 'squares', 'stars'};

makeDataNew;
FitParams = [];
plotPsycho = 1;


for run = 1:length(fittingStims)
    %% run the data set used for the psychometric function first
    % (it will run again in the next loop but whatever).
    fitting = 1;
    if strcmpi(fittingStims{run},'malania')
        
        ThreshElev = [dataStruct.(matlab.lang.makeValidName('malaniaEqual'))(1,:), ...
            dataStruct.(matlab.lang.makeValidName('malaniaLong'))(1,:), dataStruct.(matlab.lang.makeValidName('malaniaShort'))(1,:)];
        outputs = zeros(1,6);
        outputs(1:2) = aaronModel('malaniaEqual');
        outputs(3:4) = aaronModel('malaniaLong');
        outputs(5:6) = aaronModel('malaniaShort');
        FitParams = FitSigmoidArbMaxMin([outputs 1],[ThreshElev 11.5],100*ones(1,length(ThreshElev(1,:))+1)); % i added the vernier alone data here.
        
        if plotPsycho
            figure(1000)
            x = min(outputs):0.005:1;
            s=SigmoidML(FitParams,x);
            plot(x, s)
            hold on
            scatter([outputs 1], [ThreshElev 11.5]);
            hold off
            title(['Psychometric function for ', fittingStims{run}])
            saveas(gcf,['results/psychFunction,', fittingStims{run}, '.jpg'])
        end
        
    else
        ThreshElev = dataStruct.(matlab.lang.makeValidName(fittingStims{run}))(1,:);
        outputs = aaronModel(fittingStims{run});
        FitParams = FitSigmoidArbMaxMin([outputs 1],[ThreshElev 150],100*ones(1,length(ThreshElev(1,:))+1)); % i added the vernier alone data here.
        
        if plotPsycho
            figure(1000)
            x = min(outputs):0.005:1;
            s=SigmoidML(FitParams,x);
            plot(x, s)
            hold on
            scatter([outputs 1], [ThreshElev 150]);
            hold off
            title(['Psychometric function for ', fittingStims{run}])
            saveas(gcf,['results/psychFunction,', fittingStims{run}, '.jpg'])
        end
    end
    fitting = 0;

    %% Now run all conditions with the computed psychometric function
    for i = 1:length(stims)
        aaronModel(stims{i});
    end
end

