% Author:        Natacha Comandante Lou (natacom@umich.edu)
% Copyright 2019, All Rights Reserved
%
% For paper "Phenotype-Based Probabilistic
% Analysis of Heterogeneous Responses to Cancer Drugs and Their Combination
% Efficacy" by N.Comandante-Lou, M.Khaliq, D.Venkat, M.Manikkam,
% M.Fallahi-Sichani
% ------------------------------------------------------------------------
% Description: Main script for systematic simulations of varied P_death and P_stasis
% in Figure 1D

clear all; clc; close all;
N0 = [1000 0]; %initial condition : [Nlive (t=0), Ndeath (t=0)]
pdeath_range = [0,10.^[-4:0.5:-1]]; % range of P_death to be simulated
pstasis_range = [0:0.2:1];  % range of P_stasis to be simulated
[X,Y] = meshgrid([-log(1-pdeath_range),10],pstasis_range); % P_death = 1-exp(-kdeath*dt) (equation 1), quantified per hour
P_death = reshape(1-exp(-X),[],1);
P_stasis = reshape(Y,[],1);


simfoldername = strcat('./SimData/');
Home = pwd;
mkdir(simfoldername);


%% Simulation
clc;
simDate = date;
bt = 12; %time-bin size for input rate parameters(hr)
tspan = [0 bt*9]; %simulation time span
tedges = [tspan(1):bt:tspan(2)]; % time-bin for input rate parameters (hr)
ndt = 900; % number of time steps in the simulation output


nSim = 30; % number of simulations
for n = 1:nSim
    for d = 1:size(X,2)
        for s  = 1:size(X,1)
            lambda_death = X(s,d); %k_death (drug) (hr^-1)
            lambda_divi_ctrl = 0.025; %k_divi_(no drug) (hr^-1)
            
            p_stat = Y(s,d); %P_stasis
            
            %Note: In this analysis we assume input k_divi and k_death are time-invariant,
            %so to adapt them to our non-homogeneous poisson process (NHPP) simulation algorithm,
            %the input rates for the algorithm are (9 x 1) vectors of constant values instead of scalars.
            %The input rate vectors are organized in the struture K below.
            k_input.rates(:,1) = lambda_death*ones(length(tedges)-1,1);
            k_input.rates(:,2) = lambda_divi_ctrl*ones(length(tedges)-1,1);
            k_input.bedges = tspan(1):bt:tspan(end);
            k_input.p_stat = p_stat; % P_stasis
            
            %Input stoichiometry matrice for simulation of single cell reactions
            %  Death reaction: 1 Cell_Live --> 1 Cell_Dead  (k_death)
            %  Division reaction: 1 Cell_Live --> 2 Cell_Live  (k_division_nodrug)
            S = [1,0;1,0]; %[ Death Reaction substrates; Division Reaction substrates]
            P = [0,1;2,0]; %[ Death Reaction products; Division Reaction products]
            
            %Running the Gillepse Alogrithm
            [simO(s,d).T,simO(s,d).X,simO(s,d).Count] = ExactStoch_NHPP(S,P,k_input,N0,tspan) ;
        end
    end
    
    %Unify the time steps of all simulation output such that they are
    %uniformly spaced in time
    [simO] = UnifyTstep(simO,tspan,ndt);
    
    cd(simfoldername)
    save(strcat('Sim_',int2str(n),'.mat'),'simDate','simO');
    cd(Home)
end

%% Calculate fraction affected (fa) in terms of conventional metrics from the simulations
% Viability
clc;
time = 96; %time of measurement
fa_viability = zeros([size(X),nSim]);
for n = 1:nSim
    ti = find(simO(1,1).Tedges==time);
    cd(simfoldername)
    load(strcat('Sim_',int2str(n),'.mat'));
    cd(Home)
    
    ctrl = simO(1,1);
    for d = 1:size(X,2)
        for s  = 1:size(X,1)
            fa_viability(s,d,n) = 1 - simO(s,d).Nlive(ti)./ctrl.Nlive(ti);
        end
    end
end

mfa_viability = mean(fa_viability,3);
sfa_viability = std(fa_viability,0,3);


%%% DIP
dip = zeros([size(X),nSim]);
for n = 1:nSim
    ti = find(simO(1,1).Tedges==time);
    cd(simfoldername)
    load(strcat('Sim_',int2str(n),'.mat'));
    cd(Home)
    
    ctrl = simO(1,1);
    knet_ctrl = log(ctrl.Nlive(ti)/ctrl.Nlive(1))/time;
    for d = 1:size(X,2)
        for s  = 1:size(X,1)
            if simO(s,d).Nlive(ti)==0
                 dip(s,d,n) = (log(1e-16/simO(s,d).Nlive(1))/time)/knet_ctrl;
            else
            dip(s,d,n) = (log(simO(s,d).Nlive(ti)/simO(s,d).Nlive(1))/time)/knet_ctrl;
            end
        end
    end
end

fa_dip = (1-dip)./max(max(1-dip));
mfa_dip = mean(fa_dip,3);
sfa_dip = std(fa_dip,0,3);

%%% GR
gr = zeros([size(X),nSim]);
for n = 1:nSim
    ti = find(simO(1,1).Tedges==time);
    cd(simfoldername)
    load(strcat('Sim_',int2str(n),'.mat'));
    cd(Home)
    
    ctrl = simO(1,1);
    knet_ctrl = log2(ctrl.Nlive(ti)/ctrl.Nlive(1))/time;
    for d = 1:size(X,2)
        for s  = 1:size(X,1)
            if simO(s,d).Nlive(ti)==0
                 gr(s,d,n) = 2.^((log2(1e-16/simO(s,d).Nlive(1))/time)/knet_ctrl)-1;
            else
            gr(s,d,n) = 2.^((log2(simO(s,d).Nlive(ti)/simO(s,d).Nlive(1))/time)/knet_ctrl)-1;
            end
        end
    end
end

fa_gr = (1-gr)./2;
mfa_gr = mean(fa_gr,3);
sfa_gr = std(fa_gr,0,3);



%% Plot
clc;
c = [[1:length(P_death)]'*0.02, 0.6*ones(length(P_stasis),1), P_stasis];
figure()
%-----------------------------fa_viability--------------------------------%
subplot(3,3,1) %P_death
logx = log10(1-exp(-X)); logy = log10(mfa_viability);
logx(isinf(logx))= -5; logy(isinf(logy)) = -5;
for j = 1:size(X,1)
    plot(logx(j,:),logy(j,:),'k--','LineWidth',0.5); hold on;
end
scatter(logx(:),logy(:),100,c,'filled');hold on;
xlabel('P_{death}');
xl= unique(P_death)'; %xtl = compose('%0.0e',xl(1:2:end));
xtl = {'0','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'};
xticklabels(xtl);
xticks(unique(logx(:,[1,2:2:end,end])))
xlim([-5 0])
ylabel('fa');
yticklabels(xtl);
yticks(unique(logx(:,[1,2:2:end,end])))
ylim([-5 0])
set(gca,'FontSize',16);

subplot(3,3,2) %P_death
for j = 1:size(Y,2)
    plot(Y(:,j),mfa_viability(:,j),'k--','LineWidth',0.5); hold on;
end
scatter(Y(:),mfa_viability(:),100,c,'filled');
xlabel('P_{stasis}' );
ylabel('fa');
title('Viability');
set(gca,'FontSize',16);

subplot(3,3,3) %P_death or stasis
plot(-5:0.5:0,-5:0.5:0,'k','LineWidth',1); hold on
%P_{death U stasis} = P_{death}+P_{division}*P_{stasis}-P_{death}*P_{division}*P_{stasis}
Z = (1-exp(-X))+Y*(1-exp(-lambda_divi_ctrl))-(1-exp(-X)).*Y*(1-exp(-lambda_divi_ctrl));
logx = log10(Z); logy = log10(mfa_viability);
logx(isinf(logx))= -5; logy(isinf(logy)) = -5;
scatter(logx(:),logy(:),100,c,'filled');
xlabel('P_{death U stasis}');
xl= logx(1,[1,2:2:end,end]); %xtl = compose('%0.0e',xl(1:2:end));
xtl = {'0','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'};
xticklabels(xtl);
xticks(xl)
xlim([-5 0])
ylabel('fa');
yticklabels(xtl);
yticks(xl)
ylim([-5 0])
set(gca,'FontSize',16);
%-----------------------------fa_DIP--------------------------------------%
subplot(3,3,4)
logx = log10(1-exp(-X)); logy = log10(mfa_dip);
logx(isinf(logx))= -5; logy(isinf(logy)) = -5;
for j = 1:size(X,1)
    plot(logx(j,:),logy(j,:),'k--','LineWidth',0.5); hold on;
end
scatter(logx(:),logy(:),100,c,'filled');
xlabel('P_{death}');
xl= unique(P_death)'; %xtl = compose('%0.0e',xl(1:2:end));
xtl = {'0','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'};
xticklabels(xtl);
xticks(unique(logx(:,[1,2:2:end,end])))
xlim([-5 0])
ylabel('fa');
yticklabels(xtl);
yticks(unique(logx(:,[1,2:2:end,end])))
ylim([-5 0])
set(gca,'FontSize',16);

subplot(3,3,5)
for j = 1:size(Y,2)
    plot(Y(:,j),mfa_dip(:,j),'k--','LineWidth',0.5); hold on;
end
scatter(Y(:),mfa_dip(:),100,c,'filled');
xlabel('P_{stasis}' );
ylabel('fa');
title('DIP')
set(gca,'FontSize',16);

subplot(3,3,6)
plot(-5:0.5:0,-5:0.5:0,'k','LineWidth',1); hold on
Z = (1-exp(-X))+Y*(1-exp(-lambda_divi_ctrl))-(1-exp(-X)).*Y*(1-exp(-lambda_divi_ctrl));
logx = log10(Z); logy = log10(mfa_dip);
logx(isinf(logx))= -5; logy(isinf(logy)) = -5;
scatter(logx(:),logy(:),100,c,'filled');
xlabel('P_{death U stasis}');
xl= logx(1,[1,2:2:end,end]); %xtl = compose('%0.0e',xl(1:2:end));
xtl = {'0','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'};
xticklabels(xtl);
xticks(xl)
xlim([-5 0])
ylabel('fa');
yticklabels(xtl);
yticks(xl)
ylim([-5 0])
set(gca,'FontSize',16);
%----------------------------fa_GR----------------------------------------%
subplot(3,3,7)
logx = log10(1-exp(-X)); logy = log10(mfa_gr);
logx(isinf(logx))= -5; logy(isinf(logy)) = -5;
for j = 1:size(X,1)
    
    plot(logx(j,:),logy(j,:),'k--','LineWidth',0.5); hold on;
end
scatter(logx(:),logy(:),100,c,'filled');
xlabel('P_{death}');
xl= unique(P_death)'; %xtl = compose('%0.0e',xl(1:2:end));
xtl = {'0','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'};
xticklabels(xtl);
xticks(unique(logx(:,[1,2:2:end,end])))
xlim([-5 0])
ylabel('fa');
yticklabels(xtl);
yticks(unique(logx(:,[1,2:2:end,end])))
ylim([-5 0])
set(gca,'FontSize',16);

subplot(3,3,8)

for j = 1:size(Y,2)
    plot(Y(:,j),mfa_gr(:,j),'k--','LineWidth',0.5); hold on;
end
scatter(Y(:),mfa_gr(:),100,c,'filled');
xlabel('P_{stasis}' );
ylabel('fa');
title('GR')
set(gca,'FontSize',16);

subplot(3,3,9)
plot(-5:0.5:0,-5:0.5:0,'k','LineWidth',1); hold on
Z = (1-exp(-X))+Y*(1-exp(-lambda_divi_ctrl))-(1-exp(-X)).*Y*(1-exp(-lambda_divi_ctrl));
logx = log10(Z); logy = log10(mfa_gr);
logx(isinf(logx))= -5; logy(isinf(logy)) = -5;
scatter(logx(:),logy(:),100,c,'filled');
xlabel('P_{death U stasis}');
xl= logx(1,[1,2:2:end,end]); %xtl = compose('%0.0e',xl(1:2:end));
xtl = {'0','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'};
xticklabels(xtl);
xticks(xl)
xlim([-5 0])
ylabel('fa');
yticklabels(xtl);
yticks(xl)
ylim([-5 0])
set(gca,'FontSize',16);

% Legend
figure()
P_death(P_death==0) = 1e-5;
scatter(log10(P_death), P_stasis,90,c,'filled');
xticklabels({'0','10^{-2.5}','1'});
xticks([-5 -2.5 0])
xlim([-5 0])
xlabel('P_{death}');
ylabel('P_{stasis}');
set(gca,'FontSize',18);

