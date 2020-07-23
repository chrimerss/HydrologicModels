function [OV, ET, Hend, Cend] = Pdm(Hpar, Bpar, Hbeg, PP, PET);
%% function [OV, ET, Hend, Cend] = Pdm(Hpar, Bpar, Hbeg, PP, PET)
%% Code to run Probability Distribution Model for 1 time step
%% (Finite Non-Leaky tank - soft threshold) 
%% 9/11/2005 Hoshin V. Gupta
%% INPUTS
%%   Hpar = Parameter --> Maximum height of soil zone tank
%%   Bpar = Parameter --> B exponent
%%   Hbeg = Initial height of water in tank
%%   PP   = Precipitation for the time-step
%%   PET  = Potential ET for the time-step
%% OUTPUTS
%%   OV   = Model outflow (Precipitation Excess)
%%   ET   = Model ET (Actual Evapotranspiration)
%%   Hend = Final height of water in tank
%%   Cend = Final storage contents in tank     
%%       Cend = Cpar*(1-(1-(Hend/Hpar))^(1+b));
%%   where b = log(1-(Bpar/2))/log(0.5)
%%--------------------------------------------------------------------------

    OV   = []; % Initialize outflow
    ET   = []; % Initialize evap
    if Bpar < 0; error('B parameter < 0'); end;
    if Bpar > 2; error('B parameter > 2'); end;
    if Bpar == 2;
        b = 10^6;
    else;
        b    = log(1-(Bpar/2))/log(0.5); % Convert from scaled B (0-2) to unscaled b (0 - Inf)
    end;
    Cpar = Hpar/(1+b); % Compute maximum capacity of soil zone
    
%--(1)--Execute model for one time step
    Cbeg = Cpar*(1-power((1-(Hbeg/Hpar)),(1+b)));     % Contents at begining
    OV2 = max(0,PP+Hbeg-Hpar);                        % Compute OV2 if enough PP
    PPinf = PP-OV2;                                   % PP that does not go to OV2
    Hint = min(Hpar,(PPinf+Hbeg));                    % Intermediate height
    Cint = Cpar*(1-power((1-(Hint/Hpar)),(1+b)));     % Intermediate contents
    OV1 = max(0, PPinf + Cbeg - Cint);                % Compute OV1
    OV = OV1 + OV2;                                   % Compute total OV
    ET = min(PET, Cint);                              % Compute ET
    Cend = Cint - ET;                                 % Final contents
    Hend = Hpar*( 1 - power((1-Cend/Cpar),1/(1+b)) ); % Final height corresponding to SMA contents
    
% End of function Pdm   
    
    