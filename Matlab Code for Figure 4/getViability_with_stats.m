function [DRmat,DR_stats,Lnorm_Nlive,Lnorm_Nlive_stats] = getViability_with_stats(SimFolder,simID,time,dose,file_prefix)
% Author:        Natacha Comandante Lou (natacom@umich.edu)
% Copyright 2019, All Rights Reserved
%
% For paper "Phenotype-Based Probabilistic
% Analysis of Heterogeneous Responses to Cancer Drugs and Their Combination
% Efficacy" by N.Comandante-Lou, M.Khaliq, D.Venkat, M.Manikkam,
% M.Fallahi-Sichani
% ------------------------------------------------------------------------
% * Description: This function estimates the growth-rate (GR) inhibition
% (developed by Hafner et al. 2016) from simulations and calculates their 
% statistics (mean and confidence interval)
%
% * INPUT:
%   - SimFolder: simulation output folder
%   - simID: simulation index
%   - time: time point to measure GR (hr)
%   - dose: dose for the condition simulated
%   - file_prefix: prefix for the simulation output file
%   
% * OUTPUT:
% - DRmat(c).mat(s,d,t): viability estimated from the s th 
%   simultion of condition c, evaluated at specified dose d and timebin t
% - DR_stats(c).mat(d,:,t): statistics of viability 
%   estimated cross all simulations in GRmat(c).mat(s,d,t). columns: mean,low CI, high CI


Home  = pwd;
CI = [0.025  0.975]; %confidence interval

% Normalized Viability Dose Response (Time-Slices)
for n = simID
    cd(SimFolder);
    load(strcat(file_prefix,int2str(n),'.mat'));
    cd(Home)
    for j = 1:size(simO,2)
        [~,dose_id] = find(ismember([simO(:,j).Dose],dose)==1);
        [~,DR] = getViability(simO(:,j),time,dose);
        DRmat(j).mat(n,:,:) = DR';
        
        nlive = reshape([simO(dose_id,j).Nlive],[],length(dose_id));
        Lnorm_Nlive(j).mat(n,:,:) = log2(nlive./nlive(1,:)); %log2 normalized (to t=1) Nlive;
        
    end
    
end

for j = 1:length(DRmat)
    for t = 1:length(time)
        DR_stats(j).mat(:,1,t) = nanmean(DRmat(j).mat(:,:,t),1);
        DR_stats(j).mat(:,2:3,t) = getCI(DRmat(j).mat(:,:,t),CI);
        
    end
    for d = 1:length(dose_id)
        Lnorm_Nlive_stats(j).mat(:,1,d) = nanmean(Lnorm_Nlive(j).mat(:,:,d),1);
        Lnorm_Nlive_stats(j).mat(:,2:3,d) = getCI(Lnorm_Nlive(j).mat(:,:,d),CI);
    end
    
end


end


function [X,X_norm] = getViability(O,time,dose)

X = zeros(length(time),length(dose));

[~,dose_id] = find(ismember([O.Dose],dose)==1);

n = 0;
for d = dose_id
    n = n+1;
    for t = 1:length(time)
    
            [~,ti] = min(abs(O(d).Tedges-time(t))); %index of "O.T" where it is closest to the timepoint specified in "time";
            X(t,n) = O(d).Nlive(ti);
        
    end
end

X_norm = X./X(:,1);
end

function [x_CI] = getCI(x,CI)
SEM = nanstd(x,0,1)./sqrt(size(x,1)); % Standard Error
ts = tinv(CI,size(x,1)-1);          % T-Score
x_CI(:,1) = nanmean(x,1) + ts(1).*SEM;      % Lower Confidence Intervals
x_CI(:,2) = nanmean(x,1) + ts(2).*SEM;      % Upper Confidence Intervals
end