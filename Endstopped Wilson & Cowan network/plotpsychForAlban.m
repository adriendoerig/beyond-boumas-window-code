% Plots the psychometric function and the data. May need some fiddling
% depending on the data set.

modelOutputs = [1 2 3];
experimentalThresholds = [1.3 2 2.8];
fitParams = FitSigmoidArbMaxMin(modelOutputs,experimentalThresholds,100*ones(1,length(ThreshElev(1,:))));

x = 0:0.005:2;
s=SigmoidML(fitParams,x);
plot(x, s)
hold on
scatter(CrossCorr, Data(1,2:end));
scatter(CrossCorr(1), Data(1,2), 'filled');
legend('psych fct from octagons', 'whithout grouping flankers', 'with grouping flankers');