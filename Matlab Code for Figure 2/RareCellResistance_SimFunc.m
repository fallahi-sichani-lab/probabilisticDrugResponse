% Author:        Natacha Comandante Lou (natacom@umich.edu)
% Copyright 2019, All Rights Reserved
%
% For paper "Phenotype-Based Probabilistic
% Analysis of Heterogeneous Responses to Cancer Drugs and Their Combination
% Efficacy" by N.Comandante-Lou, M.Khaliq, D.Venkat, M.Manikkam,
% M.Fallahi-Sichani
% ------------------------------------------------------------------------
% *Description: Simulation function low-frequency heterogeneity
% in Figure 2 B-C.
% *INPUT:
%   - Case_ID: simulation condition id
%   - Home: home directory
%   - simI: input structure containing parameters for simulated conditions
%   - nSim: number of simulations for each condition
% *OUTPUT:
%   - O: simulation output structure, organized as follows:
%        rows: dose/treatment conditions;
%        columns: sensitive subpopulation, resistant subpopulation, homogeneously sensitive
%                 population,total population (i.e. sum of sensitive and resistant subpopulation)


function [O,K] = RareCellResistance_SimFunc(Case_ID,Home,simI,nSim)
% Folder for stochastic simulation output for each condition
simfoldername = strcat('./SimData/','Simulated_Data_Case_',int2str(Case_ID));
mkdir(simfoldername);
s = what(simfoldername);
SimDataFolder = s.path;
% Folder for figures
figfoldername = './Figures';
mkdir(figfoldername);
s = what(figfoldername);
FigureFolder = s.path;
% Folder for output metrics calculated from the simulation output
metricsfoldername = strcat('./SimOutputData/','OutputMetrics_Case_',int2str(Case_ID));
mkdir(metricsfoldername);
s = what(metricsfoldername);
SimOutputDataFolder = s.path;

file_prefix = 'Rare_Resistance_Sim_';

% Passing input parameters from the input structure, simI
dose = simI.Dose; % arbitrary, just representing with or without treatment
bt = simI.bt; % time-bin size for input rate parameters(hr)
nbin = simI.nbin; % number input time bins
tspan = [0 bt*nbin];
tedges = [tspan(1):bt:tspan(2)]; % time-bin edges for the input rates
ndt = simI.ndt;  %number of output timesteps
w = simI.w; %initial fraction of resistant population


Xic_total = simI.Xic_total; % initial cell number
Xic(1,:) = [floor(Xic_total*(1-w)),0]; %initial sensitive population: [N_live, N_death]
Xic(2,:) = [floor(Xic_total*w),0]; %initial resistant population: [N_live, N_death]

%Create the rate vectors based on the specified input values.

%Note: In this analysis we assume input k_divi and k_death are time-invariant,
%so to adapt them to our non-homogeneous poisson process (NHPP) simulation algorithm,
%the input rates for the algorithm are (nbin x 1) vectors of constant values instead of scalars.
%The input rate vectors are organized in the struture K below.

%-------------Cell type 1 (sensitive) Parameters-------------%
K(1).lamda_death(:,1) = simI.lamda_death(1).ctrl*ones(length(tedges)-1,1); % k_(S)death(no drug)
K(1).lamda_death(:,2) = simI.lamda_death(1).drug*ones(length(tedges)-1,1); % k_(S)death(drug)
K(1).lamda_divi_ctrl = repmat(simI.lamda_divi(1).ctrl*ones(length(tedges)-1,1),1,length(dose)); % k_(S)division(no drug)
K(1).p_stat = [simI.p_stasis(1).ctrl, simI.p_stasis(1).drug]; %[P_(S)stasis(no drug),P_(S)stasis(drug)]
%-------------Cell type 2 (resistant) Parameters--------------%
K(2).lamda_death(:,1) = simI.lamda_death(2).ctrl*ones(length(tedges)-1,1); % k_(R)death(no drug)
K(2).lamda_death(:,2) = simI.lamda_death(2).drug*ones(length(tedges)-1,1); % k_(R)death(drug)
K(2).lamda_divi_ctrl = repmat(simI.lamda_divi(2).ctrl*ones(length(tedges)-1,1),1,length(dose));% k_(S)division(no drug)
K(2).p_stat = [simI.p_stasis(2).ctrl, simI.p_stasis(2).drug];%[P_(S)stasis(no drug),P_(S)stasis(drug)]

ncelltype = length(K); %number of cell type
condition = {'Sensitive Alone','Resistant Alone','Homogeneous Sensitive Response','Heterogeneous Population'}; % names of the conditions to be simulated

%% Simulation
clc;

for n = 1:nSim
    disp(sprintf('Sim %d',n));
    clear simO
    for j = 1:ncelltype + 1
        
        for d = 1:length(dose)
            if j<= ncelltype % simulate each subpopulation separately
                % k_input.rates - Columns: Death Rate, Division Rate (no drug) (hr^-1); Rows: "nbin" time bins of size "bt" (hr)
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
                [simO(d,j).T,simO(d,j).X,simO(d,j).Count] = ExactStoch_NHPP(S,P,k_input,Xic(j,:),tspan) ;
                toc
                
                %label the condition in the output structure, simO
                simO(d,j).Dose = dose(d);
                simO(d,j).condition = condition{j};
                
            elseif j > ncelltype % simulate homogeneously sensitive response
                
                k_input.rates(:,1) = K(1).lamda_death(:,d); 
                k_input.rates(:,2) = K(1).lamda_divi_ctrl(:,d);
                k_input.bedges = tspan(1):bt:tspan(end);
                k_input.p_stat = K(1).p_stat(d);
                
                S = [1,0;1,0]; %[ Death Rxn substrates; Division Rxn substrates]
                P = [0,1;2,0]; %[ Death Rxn products; Division Rxn products]
                
                %Running the Gillepse Alogrithm
                tic
                [simO(d,j).T,simO(d,j).X,simO(d,j).Count] = ExactStoch_NHPP(S,P,k_input,sum(Xic,1),tspan) ;
                toc
                
                %label the condition in the output structure, simO
                simO(d,j).Dose = dose(d);
                simO(d,j).condition = condition{j};
            end
        end
        
    end
    %Unify the time steps of all simulation output such that they are
    %uniformly spaced in time
    [simO] = unifytsteps(simO,tspan,ndt);
    
    %Sum the events from the sensitive (simO(:,1))and resistant (simO(:,2))subpopulation to
    %create the heterogeneous population responses. Output for this is
    %stored in simO(:,4)
    col = size(simO,2);%number of simulation conditions
    for d = 1:length(dose)
        
        simO(d,col+1).Nlive = 0;
        simO(d,col+1).dndivi = 0;
        simO(d,col+1).dndeath = 0;
        for j = 1:ncelltype
            simO(d,col+1).Tedges = simO(d,j).Tedges;
            simO(d,col+1).Dose = simO(d,j).Dose;
            simO(d,col+1).Nlive = simO(d,col+1).Nlive + simO(d,j).Nlive;
            
            simO(d,col+1).dndivi = simO(d,col+1).dndivi + simO(d,j).dndivi;
            simO(d,col+1).dndeath = simO(d,col+1).dndeath + simO(d,j).dndeath;
        end
        
        simO(d,col+1).condition = condition{col+1};
    end
    
    simI.Xic = Xic;
    simI.Dose = dose;
    simDate = date;
    
    %Save the simulation input and output for all conditions in 'SimData'
    cd(SimDataFolder)
    save(strcat(file_prefix,int2str(n),'.mat'),'simDate','simO','simI','K','Case_ID');
    cd(Home)
    
end

%% Calculate different metrics from the simulation output

% Probabilistic Phenotype Metrics (k's)
clc;
simID = 1:nSim;
time = simO(1).Tedges;
dt = mean(diff(time));
dbin = 12; %(hr) bin size for k's estimation

[O.kmat,O.kmat_stats,~,~,O.tbins] = getKRates_with_stats(Home,SimDataFolder,simID,dose,dbin,dt,tspan,file_prefix);

% GR
[O.GRmat,O.GR_stats] = getGR_with_stats(SimDataFolder,simID,O.tbins,dose,file_prefix);

% DIP
[O.DIPmat,O.DIP_stats] = getDIP_with_stats(SimDataFolder,simID,O.tbins,dose,file_prefix);
% Viability
[O.DRmat,O.DR_stats,O.Lnorm_Nlive,O.Lnorm_Nlive_stats] = getViability_with_stats(SimDataFolder,simID,O.tbins,dose,file_prefix);

%% Save Output Metrics
O.time = time;
O.resistant_fold = K(1).p_stat(2)/K(2).p_stat(2);
O.resistant_fraction = w;

cd(SimOutputDataFolder)
save(strcat(file_prefix,'OutputMetrics_Case',int2str(Case_ID),'.mat'),'simDate','O','simI','K','Case_ID');

cd(Home)

end

%% Function called from above
function [O] = unifytsteps(O,tspan,ndt)
%This function converts event times (output from the gillepie algorithm)
%into the number of events occured in each uniformly sized small interval dt
%(dt=tspan/ndt)

[row,col] = size(O);
for j = 1:row
    for k = 1:col
        Tedges = linspace(tspan(1),tspan(2),ndt+1);
        
        dt = mean(diff(Tedges));
        Tedges(1) = Tedges(1) + 1e-5*dt;
        
        [O(j,k).dndivi,~] = histcounts(O(j,k).Count.tdivi(2:end),Tedges);
        [O(j,k).dndeath,~] = histcounts(O(j,k).Count.tdeath(2:end),Tedges);
        Xic = O(j,k).X(1,1);
        O(j,k).Nlive = Xic + cumsum(O(j,k).dndivi - O(j,k).dndeath);
        O(j,k).Nlive = [Xic, O(j,k).Nlive];
        O(j,k).Tedges = Tedges;
        
    end
end
end

