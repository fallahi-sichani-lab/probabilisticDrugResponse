% Author:        Natacha Comandante Lou (natacom@umich.edu)
% Copyright 2019, All Rights Reserved
%
% For paper "Phenotype-Based Probabilistic
% Analysis of Heterogeneous Responses to Cancer Drugs and Their Combination
% Efficacy" by N.Comandante-Lou, M.Khaliq, D.Venkat, M.Manikkam,
% M.Fallahi-Sichani
% ------------------------------------------------------------------------
% Description: Main script for systematic simulations of drug combination
% interactions in Figure 4C

clear all; close all; clc;

MetaSim_ID = 1;
Home = pwd;
SimDataFolder = strcat('MetaSim_Data_Case_',int2str(MetaSim_ID));
FigureFolder = strcat('Figures');
mkdir(strcat('./',SimDataFolder));
mkdir(strcat('./',FigureFolder));

bt = 12; % time-bin size for input rate parameters(hr)
tspan = [0 120]; %duration of the simulation
tedges = [tspan(1):bt:tspan(2)]; % time-bin for input rate parameters (hr)
ndt = 900; % number of time steps in the simulation output
Xic = [2000,0]; %initial cell number


dose = 10.^[-6, -1]; % arbitrary, just representing with or without treatment

%Create the rate vectors based on the specified input values.

%Note: In this analysis we assume input k_divi and k_death are time-invariant,
%so to adapt them to our non-homogeneous poisson process (NHPP) simulation algorithm,
%the input rates for the algorithm are (nbin x 1) vectors of constant values instead of scalars.
%The input rate vectors are organized in the struture K below.
%--------------Drug A Parameters--------------%
K(1).lamda_death(:,1) = 0*ones(length(tedges)-1,1); % k_death(no drug) = 0 (hr^-1) , input death rate of the untreated control
K(1).lamda_death(:,2) = 0.01*ones(length(tedges)-1,1); % k_death(drug A) = 0.01 (hr^-1) , input death rate induced by drug A
K(1).lamda_divi_ctrl = repmat(0.035*ones(length(tedges)-1,1),1,length(dose));% k_division(no drug) = 0.035 (hr^-1) , input division rate of the untreated control
K(1).p_stat = [0 0.2];%[P_stasis(no drug), P_stasis(drug A)]
%--------------Drug B Parameters--------------%
K(2).lamda_death(:,1) = 0*ones(length(tedges)-1,1); % k_death(no drug) = 0 (hr^-1) , input death rate of the untreated control
K(2).lamda_death(:,2) = 0.01*ones(length(tedges)-1,1); % k_death(drug B) = 0.01 (hr^-1) , input death rate induced by drug B
K(2).lamda_divi_ctrl = repmat(0.035*ones(length(tedges)-1,1),1,length(dose)); % k_division(no drug) = 0.035 (hr^-1) , input division rate of the untreated control
K(2).p_stat = [0 0.2];%[P_stasis(no drug), P_stasis(drug B)]

%Simulation Conditions
%Q.death(Q.stasis) are the model parameters that describe the degrees of
%cytotoxic (cytostatic) interactions between drug A and B relative to
%statistical independence. Mathematically, Q is the inverse of the
%combination index (CI) defined in the manuscript (Equation 12, 13).


% Q.death = [0,2.^[-2:0.5:2]]; %drug combination interaction in death
% Q.statsis = [0,2.^[-2:0.5:2]]; %drug combination interaction in stasis

Q.death = [2.^[-2:0.5:2]]; %drug combination interaction in death
Q.statsis = [2.^[-2:0.5:2]]; %drug combination interaction in stasis

nSim = 10; %number of simulations

%% Simulation
for q1 = 1:length(Q.death)
    for q2 = 1:length(Q.statsis)
        q.death = Q.death(q1);
        q.statsis = Q.statsis(q2);
        for n = 1:nSim
            disp(sprintf('Sim %d',n));
            for j = 1:length(K)+1
                for d = 1:length(dose)
                    if j <= length(K) %Single drug simulations
                        % k_input.rates - Columns: Death rate (with drug) , Division rate (no drug) (hr^-1); Rows: time bins of size bt (hr)
                        k_input.rates(:,1) = K(j).lamda_death(:,d);
                        k_input.rates(:,2) = K(j).lamda_divi_ctrl(:,d);
                        k_input.bedges = tspan(1):bt:tspan(end);
                        k_input.p_stat = K(j).p_stat(d);
                        
                        %Input stoichiometry matrice for simulation of single cell reactions
                        %  Death reaction: 1 Cell_Live --> 1 Cell_Dead  (k_death)
                        %  Division reaction: 1 Cell_Live --> 2 Cell_Live  (k_division_nodrug)
                        S = [1,0;1,0]; %[ Death Reaction substrates; Division Reaction substrates]
                        P = [0,1;2,0]; %[ Death Reaction products; Division Reaction products]
                        
                        %Running the Gillepse Alogrithm
                        tic
                        [simO(d,j).T,simO(d,j).X,simO(d,j).Count] = ExactStoch_NHPP(S,P,k_input,Xic,tspan) ;
                        toc
                        
                        %Label the condition in the output structure, simO
                        simO(d,j).Dose = dose(d);
                    end
                    
                    if j == length(K)+1 %Drug A+ Drug B Simulation
                        % k_input.rates - Columns: Death rate (induced by drug A), Death rate (induced by drug B), Division rate (no drug) (hr^-1); Rows: time bins of size bt (hr)
                        k_input.rates= [K(1).lamda_death(:,d) , K(2).lamda_death(:,d) , K(1).lamda_divi_ctrl(:,d)];
                        k_input.bedges = tspan(1):bt:tspan(end);
                        k_input.p_stat = [K(1).p_stat(d), K(2).p_stat(d)]; %[P_stasis(drug A), P_stasis(drug B)]
                        
                        %Input stoichiometry matrice for simulation of single cell reactions
                        %  Death reaction by drug A
                        %  1 Cell_Live --> 1 Cell_Dead  (k_death_drugA)
                        %  Death reaction by drug B
                        %  1 Cell_Live --> 1 Cell_Dead  (k_death_drugB)
                        %  Division reaction:
                        %  1 Cell_Live --> 2 Cell_Live  (k_division)
                        
                        S = [1,0;1,0;1,0]; %[ Death by drug A Reaction substrates; Death by drug B Reaction substrates; Division Reaction substrates]
                        P = [0,1;0,1;2,0]; %[ Death by drug A Reaction products; Death by drug B Reaction products; Division Reaction products]
                        
                        cd(Home)
                        tic
                        [simO(d,j).T,simO(d,j).X,simO(d,j).Count] = ExactStoch_NHPP_DrugCombo(S,P,k_input,q,Xic,tspan) ;
                        toc
                        simO(d,j).Dose = dose(d);
                        
                    end
                end
                
            end
            simI.Xic = Xic; %initial condition
            simI.Dose = dose;
            simDate = date;
            [simO] = UnifyTstep_Combo(simO,tspan,ndt);
            cd(SimDataFolder)
            save(strcat('MetaSim_',int2str(n),'.mat'),'simDate','simO','simI','K','MetaSim_ID');
            
            cd ..
        end
        
        
        %% Measure metrics from simulations
        simID = 1:nSim;
        time = [48]; % (hr)time of measurement
        file_prefix = 'MetaSim_';
        
        % GR
        
        [GRmat(q1,q2).s,GR_stats(q1,q2).s] = getGR_with_stats(SimDataFolder,simID,time,dose,file_prefix);
        GRmat(q1,q2).q = q;
        GR_stats(q1,q2).q = q;
        % Viability
        
        [DRmat(q1,q2).s,DR_stats(q1,q2).s,~,~] = getViability_with_stats(SimDataFolder,simID,time,dose,file_prefix);
        DRmat(q1,q2).q = q;
        DR_stats(q1,q2).q = q;
        % Probabilistic Phenotype Metrics (k's)
        dt = mean(diff(simO(1).Tedges));
        dbin = 12;
        
        [kmat(q1,q2).s,kmat_stat(q1,q2).s,~,~,~] = getKRates_with_stats(Home,SimDataFolder,simID,dose,dbin,dt,tspan,file_prefix);
        kmat(q1,q2).q = q;
        kmat_stat(q1,q2).q = q;
        
    end
end

cd(SimDataFolder)
save drug_combo_metasim_workspace
cd ..

%% Calculate Bliss Independence and Combination index
bin = floor(time/dbin); %selected time bin to plot for k
dd = 2; %drug treated
for q1 = 1:length(Q.death)
    for q2 = 1:length(Q.statsis)
        %mean of k's  for drug A+B with parameter Q.death(q1), Q.divi(q2), at
        %dose(dd), at time dbin*bin
        k_combo.divi(q1,q2) = squeeze(kmat_stat(q1,q2).s(3).divi(bin,1,dd));
        k_combo.death(q1,q2) = squeeze(kmat_stat(q1,q2).s(3).death(bin,1,dd));
        k_combo.divi_ctrl(q1,q2) = squeeze(kmat_stat(q1,q2).s(3).divi(bin,1,1));
        %mean of GR  for drug A+B with parameter Q.death(q1), Q.divi(q2), at
        %dose(dd) at time
        GR_combo(q1,q2) = squeeze(GR_stats(q1,q2).s(3).mat(dd,1));
        
        %mean of viability  for drug A+B with parameter Q.death(q1), Q.divi(q2), at
        %dose(dd) at time
        DR_combo(q1,q2) = squeeze(DR_stats(q1,q2).s(3).mat(dd,1));
        
        
        %Evaluate Bliss Independence using output from single-drug
        %simulations
        
        %Bliss Independence calcuated using probabilistic phenotype metrics
        k_indep.death(q1,q2) = squeeze(kmat_stat(q1,q2).s(1).death(bin,1,dd)) + squeeze(kmat_stat(q1,q2).s(2).death(bin,1,dd));
        k_indep.divi(q1,q2) = squeeze(kmat_stat(q1,q2).s(1).divi(bin,1,dd))*squeeze(kmat_stat(q1,q2).s(2).divi(bin,1,dd))/squeeze(kmat_stat(q1,q2).s(1).divi(bin,1,1));
        k_indep.p_stasis(q1,q2) = 1-(k_indep.divi(q1,q2)./squeeze(kmat_stat(q1,q2).s(1).divi(bin,1,1)));
        
        %Bliss Independence calcuated using fa_GR
        GR_e1 = (1-squeeze(GR_stats(q1,q2).s(1).mat(dd,1)))./2; %fa_GR(A)
        GR_e2 = (1-squeeze(GR_stats(q1,q2).s(2).mat(dd,1)))./2; %fa_GR(B)
        GR_bliss(q1,q2) = 1- (GR_e1 + GR_e2 - GR_e1.*GR_e2)*2;
        
        %Bliss Independence calcuated using fa_viability
        DR_bliss(q1,q2) = squeeze(DR_stats(q1,q2).s(1).mat(dd,1)).*squeeze(DR_stats(q1,q2).s(2).mat(dd,1)); 
        
        
    end
end

% Calculate Bliss combination index for fa_GR and fa_viability
for q1 = 1:length(Q.death)
    for q2 = 1:length(Q.statsis)
        
        CI_GR(q1,q2) = (1-GR_bliss(q1,q2))./(1-GR_combo(q1,q2));
        CI_DR(q1,q2) = (1-DR_bliss(q1,q2))./(1-DR_combo(q1,q2));
    end
end

%% Plot
close all; clc;

[X,Y] = meshgrid(Q.death,Q.statsis);

CI_GR(:,1)=NaN; % ignore conditions where q1 or q2 = 0;
CI_DR(:,1)=NaN;

color = [1-0.5*X(:)-0.1, 1-0.3*X(:),  0.5*Y(:)];
ic = X(:)==Y(:);
subplot(2,2,1)
scatter(log2(1./X(:)), log2(CI_DR(:)),50,color,'filled'); hold on;
scatter(log2(1./X(ic)), log2(CI_DR(ic)),80,'r','filled'); hold on;
plot([-2:0.5:2],[-2:0.5:2],'k-','LineWidth',1);
xlabel('log2[CI_{death}]'); ylabel('log2[ CI_{fa (viability)}]');title('Viability');
% xlim([-2 2]); ylim([-2 2]);
set(gca,'FontSize',14);

subplot(2,2,2)
scatter(log2(1./Y(:)), log2(CI_DR(:)),50,color,'filled'); hold on;
scatter(log2(1./Y(ic)), log2(CI_DR(ic)),80,'r','filled'); hold on;
plot([-2:0.5:2],[-2:0.5:2],'k-','LineWidth',1);
xlabel('log2[CI_{stasis}]'); ylabel('log2[ CI_{fa (viability)}]');title('Viability');
% xlim([-2 2]); ylim([-2 2]);
set(gca,'FontSize',14);

subplot(2,2,3)
scatter(log2(1./X(:)), log2(CI_GR(:)),50,color,'filled'); hold on;
scatter(log2(1./X(ic)), log2(CI_GR(ic)),80,'r','filled'); hold on;
plot([-2:0.5:2],[-2:0.5:2],'k-','LineWidth',1);
xlabel('log2[CI_{death}]'); ylabel('log2[ CI_{fa (GR)}]');title('GR');
% xlim([-2 2]); ylim([-2 2]);
set(gca,'FontSize',14);

subplot(2,2,4)
scatter(log2(1./Y(:)), log2(CI_GR(:)),50,color,'filled'); hold on;
scatter(log2(1./Y(ic)), log2(CI_GR(ic)),80,'r','filled'); hold on;
plot([-2:0.5:2],[-2:0.5:2],'k-','LineWidth',1);
xlabel('log2[CI_{stasis}]'); ylabel('log2[ CI_{fa (GR)}]');title('GR');
% xlim([-2 2]); ylim([-2 2]);
set(gca,'FontSize',14);

% Legend
figure();
scatter(log2(1./X(:)),log2(1./Y(:)),80,color,'filled');hold on;
scatter(log2(1./X(ic)), log2(1./Y(ic)),80,'r','filled'); 
xlabel('log2(CI_{death})'); ylabel('log2(CI_{stasis})');
set(gca,'FontSize',14);


