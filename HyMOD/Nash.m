function [Out, Xend] = Nash(K, N, Xbeg, Inp);
%% function [Out, Xend] = Nash(K, N, Xbeg, Inp)
%% Code to run Nash Cascade Model for 1 time step
%% (Series of Infinite Linear Leaky tanks) 
%% 9/18/2005 Hoshin V. Gupta
%% INPUTS
%%   K     = Parameter --> Leakage rate of the tanks
%%   N     = Parameter --> Number of tanks in series
%%   Xbeg  = Initial contents of water in the tanks
%%   Inp   = Input to first tank in the series
%% OUTPUTS
%%   Out   = Outflow from last tank in series
%%   Xend  = Final storage contents in tank     
%%--------------------------------------------------------------------------


%--(1)--Initialization
    Out  = []; % Initialize outflow
    
%--(2)--Computations
    for Res = 1:N;
        OO(Res) = K*Xbeg(Res);
        Xend(Res) = Xbeg(Res) - OO(Res);
        if Res==1; Xend(Res) = Xend(Res) + Inp; 
        else; Xend(Res) = Xend(Res) + OO(Res-1);
        end;
    end;
    Out = OO(N);

% End of function Nash