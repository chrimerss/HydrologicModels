function Model = Hymod01(Data, Period, Pars, InState)
%% function Model = Hymod01(Data, Period, Pars, InState)
%% Code to run the HyMod01 Model
%% 9/18/2005 Hoshin V. Gupta
%% INPUTS
%%   Data = Data Structure
%%   Period = Simulation Period array
%%   Pars = Parameter Structure
%%   InState = Initial State structure
% OUTPUTS
%%   Model = Model Computed variables structure - final values
%%=========================================================================

%--(1)--Initialize variables
% Data
    PP  = Data.PP;              % Precipitation flux data 
    PET = Data.PET;             % PET flux data
% Soil Moisture Accounting tank
    Hpar = Pars.Huz;             % Max height of soil moisture accounting tank
    Bpar = Pars.B;               % Distribution function shape parameter
    b = log((1-Bpar/2))/log(0.5);  % Transformation from bounded 'B' to unbounded 'b'
    Cpar = Hpar/(1+b);           % Maximum capacity of soil moisture accounting tank
    XCuz(1)  = InState.XCuz;     % Initial content of soil moisture accounting tank
    XHuz(1) = Hpar*( 1 - power((1-XCuz(1)/Cpar),1/(1+b)) ); % Initial height corresponding to SMA tank contents
% Routing tanks 
    Alp = Pars.Alp;             % Quick-slow split parameter
    Nq  = Pars.Nq;              % Number of quickflow routing tanks
    Kq  = Pars.Kq;              % Quickflow routing tanks rate parameter
    Ks  =  Pars.Ks;             % Slowflow routing tanks rate parameter
    Xq(1,1:Nq) = InState.Xq;    % Initial state values for Quickflow routing tanks
    Xs(1)      = InState.Xs;    % Slowflow routing tank   

%--(2)--Run Model Simulation TIME Loop
    Nday = length(Period);       
    for i = 1:Nday;
        Day = Period(i);
        % Run Pdm soil moisture accounting including evapotranspiration
        [OV(i), ET(i), XHuz(i), XCuz(i)] = Pdm(Hpar, Bpar, XHuz(i), PP(Day), PET(Day));
        % Run Nash Cascade routing of quickflow
        [Qq(i), Xq(i,1:Nq)] = Nash(Kq, Nq, Xq(i,1:Nq), Alp*OV(i));
        % Run Infinite Linear tank (Nash Cascade with N=1) routing of slowflow
        [Qs(i), Xs(i)] = Nash(Ks, 1, Xs(i), (1-Alp)*OV(i));
        if i<Nday;
            XHuz(i+1) = XHuz(i);
            Xq(i+1,:) = Xq(i,:);
            Xs(i+1)   = Xs(i);
        end;
    end;
    
 %--(3)-- Finalize variables   
    Model.XHuz = XHuz;            % Model computed upper zone soil moisture tank state contents
    Model.XCuz = XCuz;            % Model computed upper zone soil moisture tank state contents
    Model.Xq   = Xq;              % Model computed quickflow tank states contents
    Model.Xs   = Xs;              % Model computed slowflow tank state contents
    Model.ET   = ET;              % Model computed evapotranspiration flux
    Model.OV   = OV;              % Model computed precipitation excess flux
    Model.Qq   = Qq;              % Model computed quickflow flux
    Model.Qs   = Qs;              % Model computed slowflow flux
    Model.QQ   = Qq + Qs;         % Model computed total streamflow flux

% End of function Hymod01