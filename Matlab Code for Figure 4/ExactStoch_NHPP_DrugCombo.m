function [t,X,count] = ExactStoch_NHPP_DrugCombo(S,P,K,q,IC,tspan)
%Implementation of Gillespie Exact Stochastic Algorithm.
% Adapted from https://faculty.washington.edu/wendyt/DSDEexact.m by Wendy Thomas
% Modified by N. Comandante-Lou to simulate non-homogeneous poisson
% processes (NHPP) for modeling drug combination responses with variable
% degress of cytotoxic and cytostatic drug-interactions
% INPUT:
% S = stoichiometry of C substrates in R reactions
% P = stoichiometry of products
% K = vector of reaction rates
% IC = initial system states [Nlive, Ndeath];
% q = inverse of combination index (see Drug_Combination_MetaSim_Main.m)
% tspan = time span of the simulation (hr)
% OUTPUT:
% t = a vector that record the times of (any) events happened
% X = Number of live cells (Nlive, left column), and dead cells (Ndead, right
%  column) as a function of t

%set the uniform random number generator
jump = 10000;
rng(sum(jump*clock),'v5uniform') 
% note: jump should be a large number, such that at each iteration a
% different random seed is more likely to be selected.

%Initialization
t(1) = tspan(1);
maxT = tspan(2);
X = IC;
count.ndivi = 0; % cumulative count of divisions until each division time 
count.tdivi = 0; % time of division
count.ndeath(1).deaths = 0; % cumulative count of deaths induced by drug A until each tdeath1
count.tdeath(1).deaths = 0; % tdeath1, time of death induced by drug A
count.ndeath(2).deaths = 0; % cumulative count of deaths induced by drug B until each tdeath2
count.tdeath(2).deaths = 0; % tdeath2, time of death induced by drug A
count.nstat = 0; % cumulative count of stasis event until each stasis time 
count.tstat = 0; % time of stasis

count.nRC = 0;
[nrxn,nspc] = size(S); %nrxn - number of reactions; %nspc - number of species


%---------------------Gillespie Algorithm ------------------------
while t(end) <= maxT && X(end,1)>0
    
    %step 1: Calculate a's (reaction rates given system state, i.e.  number of dead cells and live cells)
    tvec = sort([t(end),K.bedges]);
    idx = find(tvec == t(end), 1, 'last' )-1;
    p_stat(1) = K.p_stat(1); %P_stasis (drug A)
    p_stat(2) = K.p_stat(2); %P_stasis (drug B)
    a = K.rates;  %[Death rate (induced by drug A), Death rate (induced by drug B), Division rate (no drug)]
    for b = 1:length(K.bedges)-1
        for r = 1:nrxn
            for c = 1:nspc
                if S(r,c) == 1
                    a(b,r) = a(b,r)*X(end,c);
                elseif S(r,c) == 2
                    a(b,r) = a(b,r)*X(end,c)*(X(end,c)-1)/2;
                elseif S(r,c) == 3
                    a(b,r) = a(b,r)*X(end,c)*(X(end, c)-1)/2*(X(end,c)-2)/3;
                end
            end
        end
    end
    
    %q.death is the parameter for fold change of cell killing induced by drug A + B 
    %realative to statistical independence:
    %k_death(A+B,indep) = k_death(A) + k_death(B) . - (Equation 10)
    %k_death(A+B) = k_death(A+B,indep)*q.death = k_death(A)*q.death + k_death(B)*q.death 

    a(idx,1:2) = a(idx,1:2)*q.death;
    
    
    a0 = sum(a(idx,:)); %get the total rate at which the system changes in any way
    if a0 == 0
        X(end+1,:) = X(end,:);
        t = [t; t(end)+1];
        %break;
    end
    % Step 2: Determing time of the next event (tau)
    p1 = rand;
    tau = (1/a0)*log(1/p1); % tau is exponentially distributed with parameter a0
    if t(end)+tau > maxT
        break;
    end
    % Step 3: Determine which event it is
    p2 = rand;
    for r = 1:nrxn
        if p2*a0 <= sum(a(idx,1:r)); break; end
    end
    %Update system state (number of dead cells and live cells)
    if r==1 % death event induced by drug A
        X(end+1,:) = X(end,:) - S(r,:) + P(r,:);
        t = [t; t(end)+tau];
        count.nRC = count.nRC + 1;
        count.ndeath(1).deaths = [count.ndeath(1).deaths; count.ndeath(1).deaths(end) + 1];
        count.tdeath(1).deaths = [count.tdeath(1).deaths; t(end)];
    elseif r==2 %death event induced by drug B
        X(end+1,:) = X(end,:) - S(r,:) + P(r,:);
        t = [t; t(end)+tau];
        count.nRC = count.nRC + 1;
        count.ndeath(2).deaths = [count.ndeath(2).deaths; count.ndeath(2).deaths(end) + 1];
        count.tdeath(2).deaths = [count.tdeath(2).deaths; t(end)];
    elseif r==3 % potential division event 
        p3 = rand;
        p_stat_nonindep = (p_stat(1)+p_stat(2)-p_stat(1)*p_stat(2))*q.statsis; %P_stasis (A+B) = P_stasis(A+B, indep)*q.stasis
        if p3 <= p_stat_nonindep  %division inhibited by the drug combination
            X(end+1,:) = X(end,:); %no change in cell number
            t = [t; t(end)+tau];
            count.nRC = count.nRC + 1;
            count.nstat = [count.nstat; count.nstat(end) + 1];
            count.tstat = [count.tstat; t(end)];
        else %divides any way in the presence of drugs
            X(end+1,:) = X(end,:) - S(r,:) + P(r,:);
            t = [t; t(end)+tau];
            count.nRC = count.nRC + 1;
            count.ndivi = [count.ndivi; count.ndivi(end) + 1];
            count.tdivi = [count.tdivi; t(end)];
        end
    end
end
end


