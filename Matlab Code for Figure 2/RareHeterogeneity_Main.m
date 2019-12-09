% Author:        Natacha Comandante Lou (natacom@umich.edu)
% Copyright 2019, All Rights Reserved
%
% For paper "Phenotype-Based Probabilistic
% Analysis of Heterogeneous Responses to Cancer Drugs and Their Combination
% Efficacy" by N.Comandante-Lou, M.Khaliq, D.Venkat, M.Manikkam,
% M.Fallahi-Sichani
% ------------------------------------------------------------------------
% Description: Main script for simulations of low-frequency heterogeneity
% in Figure 2 B-C.
% 

clear all; close all; clc;
Home = pwd;

%% Simulation Conditions for Fig 2B
% R = 2.^[0 4]; % r, level of resistance
% W = [0:0.001:0.05]; % w, initial fraction of resistant population

R = 2.^[0:0.1:4];
W = [0 0.03];


for r = 1:length(R)
    for w = 1:length(W)
        Case_ID = (r-1)*length(W)+ w;
        
        simI.bt = 12; % time-bin size for input rate parameters(hr)
        simI.nbin = 10; % number input time bins
        simI.ndt = 900; %number of output timesteps
        simI.Xic_total = 300; % initial cell number 
        simI.w = W(w); %initial fraction of resistant population
        
        simI.Dose = 10.^[-6, -1]; % arbitrary, just representing with or without treatment
        
        %Input parameters for the sensitive population (S)
        simI.lamda_divi(1).ctrl = 0.035; %k_(S)division(no drug), input division rate of the untreated control (hr^-1) 
        simI.lamda_death(1).ctrl = 0; %k_(S)division(drug),input death rate of the untreated control (hr^-1)
        simI.lamda_death(1).drug = 0.03; %k_(S)death(drug), input death rate of the drug treated condition (hr^-1)  
        simI.p_stasis(1).ctrl = 0; %P_(S)stasis(no drug), input P_stasis of the untreated control
        simI.p_stasis(1).drug = 0.8;%P_(S)stasis(drug),input P_stasis of the drug treated condition
        
        %Input parameters for the resistant population (R)
        simI.lamda_divi(2).ctrl = simI.lamda_divi(1).ctrl*0.57; %k_(R)division(no drug), input division rate of the untreated control (hr^-1) 
        simI.lamda_death(2).ctrl = simI.lamda_death(1).ctrl; %k_(R)division(drug),input death rate of the untreated control (hr^-1)
        simI.lamda_death(2).drug = simI.lamda_death(1).drug/R(r); %k_(R)death(drug), input death rate of the drug treated condition (hr^-1)  
        simI.p_stasis(2).ctrl = simI.p_stasis(1).ctrl; %P_(R)stasis(no drug), input P_stasis of the untreated control
        simI.p_stasis(2).drug = simI.p_stasis(1).drug/R(r); %P_(R)stasis(drug),input P_stasis of the drug treated condition
        
        nSim = 50; %Number of simulations
        
        [O(r,w),K(r,w).k] = RareCellResistance_SimFunc(Case_ID,Home,simI,nSim); %Run Simulations

    end
end

save metasim_b

%% Simulation Conditions for Fig 2C
R = 2.^[0:0.1:4];
W = [0 0.03];


for r = 1:length(R)
    for w = 1:length(W)
        Case_ID = (r-1)*length(W)+ w;
        
        simI.bt = 12; % time-bin size for input rate parameters(hr)
        simI.nbin = 10; % number input time bins
        simI.ndt = 900; %number of output timesteps
        simI.Xic_total = 300; % initial cell number 
        simI.w = W(w); %initial fraction of resistant population
        
        simI.Dose = 10.^[-6, -1]; % arbitrary, just representing with or without treatment
        
        %Input parameters for the sensitive population (S)
        simI.lamda_divi(1).ctrl = 0.035; %k_(S)division(no drug), input division rate of the untreated control (hr^-1) 
        simI.lamda_death(1).ctrl = 0; %k_(S)division(drug),input death rate of the untreated control (hr^-1)
        simI.lamda_death(1).drug = 0.03; %k_(S)death(drug), input death rate of the drug treated condition (hr^-1)  
        simI.p_stasis(1).ctrl = 0; %P_(S)stasis(no drug), input P_stasis of the untreated control
        simI.p_stasis(1).drug = 0.8;%P_(S)stasis(drug),input P_stasis of the drug treated condition
        
        %Input parameters for the resistant population (R)
        simI.lamda_divi(2).ctrl = simI.lamda_divi(1).ctrl*0.57; %k_(R)division(no drug), input division rate of the untreated control (hr^-1) 
        simI.lamda_death(2).ctrl = simI.lamda_death(1).ctrl; %k_(R)division(drug),input death rate of the untreated control (hr^-1)
        simI.lamda_death(2).drug = simI.lamda_death(1).drug/R(r); %k_(R)death(drug), input death rate of the drug treated condition (hr^-1)  
        simI.p_stasis(2).ctrl = simI.p_stasis(1).ctrl; %P_(R)stasis(no drug), input P_stasis of the untreated control
        simI.p_stasis(2).drug = simI.p_stasis(1).drug/R(r); %P_(R)stasis(drug),input P_stasis of the drug treated condition
        
        nSim = 50; %Number of simulations
        
        [O(r,w),K(r,w).k] = RareCellResistance_SimFunc(Case_ID,Home,simI,nSim); %Run Simulations

    end
end

save metasim_c