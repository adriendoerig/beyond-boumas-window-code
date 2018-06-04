% Plots the psychometric function and the data. May need some fiddling
% depending on the data set.
x = 0:0.005:2;
s=SigmoidML(FitParams,x);
plot(x, s)
hold on
scatter(CrossCorr, Data(1,2:end));
scatter(CrossCorr(1), Data(1,2), 'filled');
legend('psych fct from octagons', 'whithout grouping flankers', 'with grouping flankers');