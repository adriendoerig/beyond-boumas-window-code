function params = FitSigmoidArbMaxMin(x,pc,n,varargin)

% FITSIGMOIDARBMAXMIN Fit a sigmoidal function to data using maximum likelihood.
% The fit sigmoidal function may have an arbitrary maximum and minimum.
% 
% Usage: [params] = FitSigmoidArbMaxMin(x,pc,n);
% 
% e.g.
% 
% x = 0:1:100;
% TrueParams = [50 -0.25 0.1 0.85];
% pc = SigmoidML(TrueParams,x);
% nz = randn(size(x))*0.05;
% nzPC = pc+nz;
% nzPC(nzPC>TrueParams(4)) = TrueParams(4);
% nzPC(nzPC<TrueParams(3)) = TrueParams(3);
% n = ones(size(x))*100;
% 
% [FitParams] = FitSigmoidArbMaxMin(x,pc,n);
% 
% X = linspace(min(x),max(x),100);
% figure('Color','w');
% plot(x,nzPC,'k.',X,SigmoidML(FitParams,X),'r-',X,SigmoidML(TrueParams,X),'b--');
% legend('Noisy Data','Fit','Model');
% 
% Created by Aaron Clarke, NDSU, Visual Dynamics Lab, February 10, 2008

% Estimate the initial parameters -------

if nargin>3
    InitParams = varargin{1};
else

    % Max and min
    TheMax = max(pc);
    TheMin = min(pc);
    
    % Mean
    dPC = abs(diff(pc));
    dX = diff(x);
    BetweenX = x(1:end-1)+dX/2;
    [MaxDiff,MeanInd] = max(dPC);
    m = BetweenX(MeanInd);
    
    % Standard deviation
    s = (m-x)./(log((TheMax-TheMin)./(pc-TheMin))-1);
    S = s(s~=0);
    nS = length(S);
    nX = length(x);
    PredMat = zeros(nS,nX);
    for iS = 1:nS
        PredMat(iS,:) = SigmoidML([m 1/S(iS) TheMin TheMax],x(:)');
    end

    Resids = sum((PredMat - repmat(pc(:)',[nS 1])).^2,2);
    [LowestResid,SInd] = min(Resids);
    s = S(SInd);
    
    if isempty(s)
        s = 1;
    end
    InitParams = [m 1/s TheMin TheMax];

end

 % Fit the data
if any(pc)
    [params,FVal,ExitFlag] = fminsearch(@(myCoef) SigmoidArbMaxMinMinFun(myCoef, x, pc, n), InitParams);
    while ExitFlag==-1  
        [params,FVal,ExitFlag] = fminsearch(@(myCoef) SigmoidArbMaxMinMinFun(myCoef, x, pc, n), params);
    end
    
else
    params = [mean(x) Inf 0 0];    
end

end




function [y] = SigmoidArbMaxMinMinFun(params, x, pc, n)

% SIGMOIDARBMAXMINMinFUN Minimization function for FitSigmoidML.m
% 
% Usage: [y] = SigmoidArbMaxMinMinFun(params, x, pc, n);
% 
% e.g.
% 
% params = [50 1/4];
% x = 0:10:100;
% n = ones(size(x))*100;
% pc = sigmoid(params,x);
% 
% y = SigmoidArbMaxMinMinFun(params, x, pc, n);
% 
% Created by Aaron Clarke, NDSU, Visual Dynamics Lab, August 25, 2008

n = n(:);
pc = pc(:);
x = x(:);

if params(3)>=min(pc) && params(4)<=max(pc) && ...
        all(SigmoidML(params, x)) && all(1-SigmoidML(params, x)) 
    y = -sum( n .* pc .*log(SigmoidML(params, x)) + ...
            n .* (1-pc) .*log(1-SigmoidML(params, x)) );       
else
    y = Inf;
end

end


