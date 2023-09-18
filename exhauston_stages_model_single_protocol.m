%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %   
% T cell exhaustion model: dynamics of a single treatment protocol      %
% Authors: Irina Kareva and Jana Gevertz                                %
% Last update: 7/7/2023                                                 %
% - Model assumes T cells can be reversibly or irreversibly exhausted,  %
%   driven by inflammation triggered by tumor killing. Drug only works  %
%   on reversibly exhausted T cells.                                    %
% - Predicts and plots dynamics of tumor, drug, T cells (in various     %
%   states), and inflammatory marker                                    %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clearvars; clc; close all;  

%% Set baseline parameters 
[p, ICs] = set_parameters_ICs(); 
tStart1 = 100; % drug 1, time of start
tmax=4000; % end time for solving ODEs

%% Set protocol
p.dose1 = 6;   % dose fixed 
p.doseD1 = p.dose1*10^3/p.V1_1; % dose recalculations fron to ng/ml
p.dosenumD1 = 8;   % number of doses fixed
p.freq1 = 8; % spacing between doses (days)
p.intvlD1 = 24*p.freq1; % dosing interval (hours)
[tDose,cDose,ti] = set_protocol(p,tmax,tStart1);

%% Solve model and process data
sol = solve_model(tDose,cDose,ICs,tmax,p,ti);
time = []; D1_p = []; x1 = []; 
y1 = []; y2 = []; y3 = []; F = [];
for k = 1:size(sol,2)
    time = [time sol{k}.x];
    D1_p = [D1_p sol{k}.y(2,:)]; 
    x1 = [x1 sol{k}.y(4,:)]; 
    y1 = [y1 sol{k}.y(5,:)]; 
    y2 = [y2 sol{k}.y(6,:)]; 
    y3 = [y3 sol{k}.y(7,:)]; 
    F = [F sol{k}.y(8,:)]; 
end


%% Plots
figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.45, 0.95]);
subplot(2,3,1,'LineWidth',2) % drug concentration
semilogy(time,D1_p,'LineWidth',2)
xlabel('Time (hours)','FontSize',14)
ylabel('Plasma drug concentration','FontSize',14) 
 
subplot(2,3,2) % tumor
plot(time,x1,'LineWidth',2)
xlabel('Time (hours)','FontSize',14)
ylabel('Tumor volume','FontSize',14)

subplot(2,3,3) % inflammation
plot(time,F,'LineWidth',2);
xlabel('Time (hours)','FontSize',14)
ylabel('Inflammation','FontSize',14)

subplot(2,3,4) % cytotoxic T cells y1
plot(time,y1,'LineWidth',2)
xlabel('Time (hours)','FontSize',14)
ylabel('Cytotoxic T cells, y_1','FontSize',14)

subplot(2,3,5) % reversibly exhausted T cells y2
plot(time,y2,'LineWidth',2)
xlabel('Time (hours)','FontSize',14)
ylabel('Reversibly exhausted T cells, y_2','FontSize',14)
      
subplot(2,3,6) % terminally exhausted T cells y3
plot(time,y3,'LineWidth',2)
xlabel('Time (hours)','FontSize',14)
ylabel('Terminally exhausted T cells, y_3','FontSize',14)

%%%%%%%%% Functions%%%%%%%%%%%%%
function [p, ICs] = set_parameters_ICs()
    %% Define variables and set initial conditions
    % 1 = D1_iv = drug 1 absorption
    % 2 = D1_p = drug 1 (pembro) in plasma
    % 3 = D1_t = drug 1 in tissue 
    % 4 = x1 = tumor
    % 5 = y1 = cytotoxic T cells
    % 6 = y2 = reversibly exhausted T cells
    % 7 = y3 = terminally exhausted T cells
    % 8 = F = inflammatory factor
    numVar = 8;
    ICs = zeros(1,numVar); % all variables initially 0 except: 
    ICs(4) = 50; % initial tumor size
    ICs(5) = 0.1; % initial number of cytotoxic T cells

    %% PK parameters
    p.V1_1  = 70; % ml/kg 
    p.V2_1  = 33;% ml/kg
    
    p.iv = 0;
    if p.iv == 1
        p.k01_1 = 0;
    elseif p.iv == 0
        p.k01_1 = 0.11; % SC
    end
    
    p.Cl1_1 = 20;% ml/kg/day
    p.Cl2_1 = 22;% ml/kg/day
    
    p.k10_1 = p.Cl1_1/p.V1_1;
    p.k12_1 = p.Cl2_1/p.V1_1;
    p.k21_1 = p.Cl2_1/p.V2_1;
  
    %% Tumor parameters
    p.lamba1 = .025;
    p.K = 3000;

    %% Immune parameters
    p.b = 0.09; %rate of cell kill by y1
    p.c = 0.01; %rate of cell kill by y3
    p.b1  = 0.1; %y1 expansion as a results of cell kill by y1
    p.xi1 = 1; % for ratio-dependent term, sclaing factor for y1 cell kill
    p.xi2 = 1; % for ratio-dependent term, sclaing factor for y3 cell kill
    p.d1  = 0.044; % normal death rate of y1

    p.g1  = 0.045; % rate of transition from y1->y2 when F>F1

    p.F1  = 10; % threshold for transition from y1->y2
    p.F2  = 25; % threshold for transition from y2->y3

    p.g2  = 0.01; % rate of return from y2 to y1 when F<F1
    p.g3  = 0.05; %rate of transition of y2->y3 when F>F2
    p.d3  = 0.01; %rate of death for terminally exhausted y3
    p.b2  = 0.02; % rate of F expansion as a function of cell kill (measure of cytokines/inflammation)
    p.d4  = 0.01; %rate of cytokine/inflammation natural clearance

    p.pd1 = 0.001; % rate of drug-mediated transition from y2 to y1
end

function [tDose,cDose,ti] = set_protocol(p,T,tStart1)
    tDose = zeros(1,p.dosenumD1); cDose = zeros(1,p.dosenumD1);
    for i=1:p.dosenumD1
        tDose(i)=p.intvlD1*(i-1) + tStart1;
        cDose(i)=p.doseD1;
    end
    % As treatment doesn't start until t = tStart = tDose(1):
    ti = [0 tDose max(T)]; 
    cDose = [0 cDose];
    tDose = [0 tDose];
end

function sol_all = solve_model(tDose,cDose,ICs,T,p,ti) 
    options_DE = odeset('RelTol',1.0e-5);
%     tPoints = []; cPoints = []; 
    for i = 1:length(tDose)
        %fprintf('Solving model for dose %d, starting at t = %d\n',i,ti(i));
        if p.iv == 1 %IV 
            if i == 1
                ICs(2) = ICs(2) + cDose(i); % D1_p0 = ICs(2);
            else
                D1_p0 = D1_p0 + cDose(i);
                ICs = [D1_iv0 D1_p0 D1_t0 x10 y10 y20 y30 F0]; 
            end
        elseif p.iv == 0 % SC
            if i == 1
                ICs(1) = ICs(1) + cDose(i); % D1_iv0 = ICs(1);
            else
                D1_iv0 = D1_iv0 + cDose(i);
                ICs = [D1_iv0 D1_p0 D1_t0 x10 y10 y20 y30 F0];
            end
        end
        tSpan = [ti(i) ti(i+1)];
        sol_all{i} = ode23s(@(t,x) ode_model(t,x,p),tSpan,ICs,options_DE);        
        D1_iv0 = sol_all{i}.y(1,end); 
        D1_p0 = sol_all{i}.y(2,end); 
        D1_t0 = sol_all{i}.y(3,end); 
        x10 = sol_all{i}.y(4,end); 
        y10 = sol_all{i}.y(5,end); 
        y20 = sol_all{i}.y(6,end);
        y30 = sol_all{i}.y(7,end); 
        F0 = sol_all{i}.y(8,end); 
    end
end


function dydt=ode_model(t,y,p)
    % Drug
    D1_iv = y(1);
    D1_p  = y(2);
    D1_t  = y(3);
    % Tumor
    x1    = y(4);
    % Immune
    y1    = y(5);
    y2    = y(6);
    y3    = y(7);
    F     = y(8);

    %% Drug 1: pembro, 2 compartment       
    dD1_iv = -p.k01_1*D1_iv;
    dD1_p = p.k01_1*D1_iv-p.k10_1*D1_p...
           -p.k12_1*D1_p+p.k21_1*D1_t*p.V2_1/p.V1_1;
    dD1_t = p.k12_1*D1_p*p.V1_1/p.V2_1-p.k21_1*D1_t;
    
    %% Tumor
    kill_term = p.b*y1/(p.xi1*x1+y1) + p.c*y3/(p.xi2*x1+y3);
    dx1 = p.lamba1*x1*(1-x1/p.K) - x1*(kill_term);

    %% Immune
    dy1 = x1*p.b1*y1/(p.xi1*x1+y1) - p.d1*y1... % growth from killing cancer cells
         -p.g1*y1*(subplus(F-p.F1)/((F-p.F1)))... % if F>F1, y1->y2
         +p.g2*y2*(subplus(p.F1-F)/((p.F1-F)))... % if F < F1, y2->y1
         +p.pd1*D1_p*y2; %rate of drug-induced inflow from y2 to y1
    dy2 = p.g1*y1*(subplus(F-p.F1)/((F-p.F1)))... % if F>F1, y1->y2
         -p.g2*y2*(subplus(p.F1-F)/((p.F1-F)))... % if F < F1, y2->y1
         -p.g3*y2*(subplus(F-p.F2)/((F-p.F2)))... %if F > F2, y2->y3
         -p.pd1*D1_p*y2; %rate of drug pushing cells from y2 to y1
    dy3 = p.g3*y2*(subplus(F-p.F2)/((F-p.F2))) - p.d3*y3;
    dF  = p.b2*x1*(kill_term)-p.d4*F; % cell-kill driven inflammation

    %%
    dydt = [dD1_iv; dD1_p; dD1_t;... % drug 1 and target 1
            dx1; dy1; dy2; dy3; dF]; % tumor and immune
end

function [tumor_end_size,tumor_erad_fails,AUC,Cmin,Cavg] = ...
    protocol_sweep(Nsweep,p,first_range,second_range,ICs,tmax,tStart1,total_drug)

    tumor_end_size = zeros(length(first_range),length(second_range)); 
    tumor_erad_fails = zeros(length(first_range),length(second_range));  
    AUC = zeros(length(first_range),length(second_range)); 
    Cmin = zeros(length(first_range),length(second_range)); 
    Cavg = zeros(length(first_range),length(second_range)); 
    for i = 1:length(first_range)
        for j=1:length(second_range)
            if Nsweep == 1
                % Set protocol
                p.dose1 = first_range(i);
                p.doseD1 = p.dose1*10^3/p.V1_1; % dose recalculations fron to ng/ml
                p.freq1 = second_range(j);
                p.intvlD1 = 24*p.freq1; % dosing interval (hours)
            elseif Nsweep == 2
                % Set protocol
                p.dosenumD1 = first_range(i);
                p.freq1 = second_range(j);
                p.intvlD1 = 24*p.freq1; % dosing interval (hours)
            elseif Nsweep == 3
                % Set protocol
                p.dosenumD1 = first_range(i);
                p.doseD1 = total_drug/p.dosenumD1;
                p.freq1 = second_range(j);
                p.intvlD1 = 24*p.freq1; % dosing interval (hours)
            end
        
            [tDose,cDose,ti] = set_protocol(p,tmax,tStart1);
            % Solve model and process data
            sol = solve_model(tDose,cDose,ICs,tmax,p,ti);
            time = []; D1_p = []; x1 = []; 
            for k = 1:size(sol,2)
                time = [time sol{k}.x];
                D1_p = [D1_p sol{k}.y(2,:)]; 
                x1 = [x1 sol{k}.y(4,:)]; 
            end
            tumor_end_size(i,j) = x1(end);
            if tumor_end_size(i,j)<1 % eradication
                tumor_erad_fails(i,j) = 1;
            else % failure
                tumor_erad_fails(i,j) = 0;
            end

            for k=1:length(time)
                if D1_p(k) < 10e-10
                    D1_p(k) = 0;
                end
            end
            [Cmin(i,j),Cavg(i,j),AUC(i,j)] = calculate_statistics(time,D1_p,ti);

%         figure; 
%         sgtitle(['Dose = ' num2str(p.dose1) ', Spacing = ' ...
%             num2str(p.intvlD1) ' hours']);
%         subplot(2,2,1) % tumor
%         plot(time,x1)
%         xlabel('Time (hours)')
%         ylabel('Tumor')
%         
%         subplot(1,2,2) % drug concentration
%         semilogy(time,D1_p)
%         xlabel('Time (hours)')
%         ylabel('Drug concentration')
            if Nsweep == 1
                fprintf('For dose = %f and frequency = %d hours:\n',p.dose1,...
                    p.intvlD1);
            elseif Nsweep == 2
                fprintf('For dose number = %d and frequency = %d hours:\n',...
                    p.dosenumD1,p.intvlD1);
            elseif Nsweep == 3
                fprintf('For dose number = %d (give %f drug per dose) and frequency = %d hours:\n',...
                    p.dosenumD1,p.doseD1*(p.V1_1/10^3),p.intvlD1);
            end
            fprintf('\tFinal tumor size = %f\n',tumor_end_size(i,j));
            fprintf('\tAUC = %f\n',AUC(i,j));
            fprintf('\tCmin = %f\n',Cmin(i,j));
            fprintf('\tCavg = %f\n',Cavg(i,j));
        end
    end
end


