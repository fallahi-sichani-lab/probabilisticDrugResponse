function [t,X,count] = ExactStoch_NHPP(S,P,K,IC,tspan)
%Implementation of Gillespie Exact Stochastic Algorithm.
% Adapted from https://faculty.washington.edu/wendyt/DSDEexact.m by Wendy Thomas
% Modified by N. Comandante-Lou to simulate non-homogeneous poisson process (NHPP)
% INPUT:
% S = stoichiometry of C substrates in R reactions
% P = stoichiometry of products
% K = vector of reaction rates
% IC = initial system states [Nlive, Ndeath];
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
count.ndivi = 0; % cumulative count of divisions until each division time (tdivi)
count.tdivi = 0;
count.ndeath = 0; % cumulative count of deaths until each death time (tdeath)
count.tdeath = 0;
count.nstat = 0; % % cumulative count of stasis event until each stasis time (tstat)
count.tstat = 0;

count.nRC = 0;
[nrxn,nspc] = size(S); %nrxn - number of reactions; %nspc - number of species

%-----------------------Gillespie Algorithm --------------------------%
while t(end) <= maxT && X(end,1)>0
    
  %Step 1: Calculate a's (reaction rates given system state, i.e.  number of dead cells and live cells)
    
    tvec = sort([t(end),K.bedges]);
    idx = find(tvec == t(end), 1, 'last' )-1; %find out the time bin where the current simulation time "t" is at 
    p_stat = K.p_stat; %P_stasis
    a = K.rates;
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
    
    a0 = sum(a(idx,:)); %get the total rate at which the system changes in any way
    
    if a0 == 0 
        X(end+1,:) = X(end,:);
        t = [t; t(end)+1];
        %break;
    end
  %Step 2: Determing time of the next event (tau)
    p1 = rand;
    tau = (1/a0)*log(1/p1); % tau is exponentially distributed with parameter a0
    if t(end)+tau > maxT
        break; % exit the simulation if the time of next event is beyond maxT
    end
  %Step 3: Determine which event it is
    p2 = rand;
    for r = 1:nrxn
        if p2*a0 <= sum(a(idx,1:r)); break; end
    end
    
    %Update system state (number of dead cells and live cells)
    if mod(r,2) == 0 % the event is a potential division event 
        p3 = rand;
        if p3 <= p_stat %being inhibited by the drug with P_stasis
            X(end+1,:) = X(end,:); %no change in cell number
            t = [t; t(end)+tau];
            count.nRC = count.nRC + 1;
            count.nstat = [count.nstat; count.nstat(end) + 1];
            count.tstat = [count.tstat; t(end)];
        else %still divides in the presence of drug with 1-P_stasis
            X(end+1,:) = X(end,:) - S(r,:) + P(r,:);
            t = [t; t(end)+tau];
            count.nRC = count.nRC + 1;
            count.ndivi = [count.ndivi; count.ndivi(end) + 1];
            count.tdivi = [count.tdivi; t(end)];
            
        end
    end
    if mod(r,2)== 1 % the event is a death event
        X(end+1,:) = X(end,:) - S(r,:) + P(r,:);
        t = [t; t(end)+tau];
        count.nRC = count.nRC + 1;
        count.ndeath = [count.ndeath; count.ndeath(end) + 1];
        count.tdeath = [count.tdeath; t(end)];
    end
    
end
end


