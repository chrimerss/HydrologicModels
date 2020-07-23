function Results(Period, Data, Model, Pars)
%% function Results(Period, Data, Model)
%% Code to Plot results after running the HyMod01 Model
%% 9/18/2005 Hoshin V. Gupta
%% INPUTS
%%   Period = Simulation Period array
%%   Data = Data Structure
%%   Model = Model Computed variables structure - initial values
% OUTPUTS
%%   None
%%=========================================================================

    DayFirst = Period(1);
    DayLast  = Period(length(Period));

    figure(1); % Open figure window number 1
    subplot(4,1,1); % Plot the Precip Data for the Period
        bar(Period, Data.PP(Period)); 
        grid on;
        axis([DayFirst DayLast 0 1.1*max(Data.PP(Period))]);
        axis ij; 
        grid on;
        ylabel('PP (mm/day)');
        title('Leaf River Daily Data'); 
        legend('P as Bar Plot'); 
    subplot(4,1,2); % Plot the Potential ET and Computed Actual ET
        plot(Period, Data.PET(Period),'r-','linewidth',1.5); 
        hold on;
        plot(Period, Model.ET,'g-','markersize',10);
        hold off;
        axis([DayFirst DayLast 0 1.1*max(max(Data.PET(Period), max(Model.ET)))]);
        grid on;
        ylabel('ET (mm/day)');
        legend('PET','AET');
    subplot(4,1,3); % Plot the Streamflow
        area(Period, Model.QQ,'edgecolor','k','facecolor',0.95*[0.6 0.9 1]);
        hold on;
        plot(Period, Data.QQ(Period),'r.','markersize',10);
        hold off; 
        axis([DayFirst DayLast 0 1.1*max(max(Data.QQ(Period),max(Model.QQ)))]);
        grid on;
        ylabel('QQ (mm/day)');
        [bias, rmse, nsce] = nanhydrostat(Data.QQ(Period),Model.QQ');
        legend(['Simulated QQ - bias(%) = ', num2str(bias), ' RMSE(%) = ', num2str(rmse), ' NSCE = ', num2str(nsce)] ,'Observed QQ');
%         legend('Simulated QQ ' ,'Observed QQ', 1);
    subplot(4,1,4); % Plot the Streamflow
        semilogy(Period, Model.QQ,'k');
        hold on;
        semilogy(Period, Data.QQ(Period),'r.','markersize',10);
        hold off; 
        axis([DayFirst DayLast min(Data.QQ(Period)) 1.1*max(max(Data.QQ(Period),max(Model.QQ)))]);
        grid on;
        ylabel('QQ (mm/day)');
        legend('Simulated QQ','Observed QQ');
        
    figure(2);
    subplot(1+Pars.Nq+1,1,1); % Plot the Soil Moisture Accounting tank state
        area(Period, Model.XCuz,'edgecolor','k','facecolor',0.95*[0.6 0.9 1]);
        axis([DayFirst DayLast 0 1.1*max(Model.XCuz)]);
        grid on;
        ylabel('XCuz (mm)');
        legend('SMA Tank state');
    for i=1:Pars.Nq
    subplot(1+Pars.Nq+1,1,i+1);  % Plot the Quickflow tanks states
        area(Period, Model.Xq(:,i),'edgecolor','k','facecolor',0.95*[0.6 0.9 1]);
        axis([DayFirst DayLast 0 1.1*max(Model.Xq(:,i))]);
        grid on;
        ylabel('Xq (mm)');
        legend('Quickflow Tank state');
    end;
    subplot(1+Pars.Nq+1,1,1+Pars.Nq+1);  % Plot the Slowflow tank state
        area(Period, Model.Xs,'edgecolor','k','facecolor',0.95*[0.6 0.9 1]);
        axis([DayFirst DayLast 0 1.1*max(Model.Xs)]);
        grid on;
        ylabel('Xs (mm)');
        legend('Slowflow Tank state');
    
    figure(3); 
        plot(Data.QQ(Period), Model.QQ,'b.','markersize',15); % Show scatterplot of obs vs sim flows
        axis('square');
        xlabel('Measured QQ (mm/day)');
        ylabel('Computed QQ (mm/day)');
        title('Scatterplot Leaf River Simulation');
        grid on;
        
% End of function Results
    
    
    