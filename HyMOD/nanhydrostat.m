function [bias, rmse, nsce] = nanhydrostat(obs, sim)
%Computes statistics commonly used in hydrologic model comparison from two
%time series "obs" and "sim": Percent Bias, Root Mean Squared Error (RMSE)
%and Nash-Sutcliff Coefficient of Efficiency (NSCE).
nancount = find(isnan(obs) == 1 | isnan(sim));
A = sim;
f = length(A) - length(nancount);
sumA = sum(A(isnan(A)==0 & isnan(obs)==0));
sumobs = sum(obs(isnan(A)==0 & isnan(obs)==0));
difAB = A(isnan(A)==0 & isnan(obs)==0) - obs(isnan(A)==0 & isnan(obs)==0);
sumdifAB = sum(difAB);
sqdifAB = difAB .^ 2;
sumE = sum(sqdifAB);
avA = mean(A(isnan(A)==0 & isnan(obs)==0));
avB = mean(obs(isnan(A)==0 & isnan(obs)==0));

%CALCULATING BIAS
%mbiasA = avA / avB;

%CALCULATING PERCENT BIAS
bias = (sumdifAB / sumobs) * 100;

%CALCULATING ROOT MEAN SQUARE ERROR PERCENTAGE
sqE = sumE / f;
RootMeSqErr = sqrt(sqE);
rmse = (RootMeSqErr / avB) * 100;

%NASH-SUTCLIFFE COEFFICIENT
difBaveB = obs(isnan(A)==0 | isnan(obs)==0) - avB;
sqdifBaveB = difBaveB .^ 2;
sumdifBaveB = sum(sqdifBaveB);
nsce = 1 - (sumE / sumdifBaveB) ;