function [GRmat,GR_stats] = getGR_with_stats(SimFolder,simID,time,dose,file_prefix)
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
% - GRmat(c).mat(s,d,t): GR metrics estimated from the s th
%   simultion of condition c, evaluated at specified dose d and timebin t
% - GR_stats(c).mat(d,:,t): statistics of GR
%   estimated cross all simulations in GRmat(c).mat(s,d,t). columns: mean,low CI, high CI

Home  = pwd;
for n = simID
    cd(SimFolder);
    load(strcat(file_prefix,int2str(n),'.mat'));
    cd(Home)
    for j = 1:size(simO,2)
        [~,GR] = getGR(simO(:,j),time,dose);
        GRmat(j).mat(n,:,:) = GR(2:end,:)';
    end
end
CI = [0.025  0.975]; %confidence interval
for j = 1:length(GRmat)
    for t = 1:length(time)
        GR_stats(j).mat(:,1,t) = nanmean(GRmat(j).mat(:,:,t),1);
        GR_stats(j).mat(:,2:3,t) = getCI(GRmat(j).mat(:,:,t),CI);
    end
end


end


function [X,GR] = getGR(O,time,dose)
% Growth Rate Inhibition (GR): normalized growth rate (at time) to untreated control
%Note:
% - GR only defined for the timepoint>0 (because log2(X_norm(1,1))=0 in the denominator);
% - Analysis is based on data with unified time steps
%INPUT:
% - O: simulation output structure
% - time: A vector that specifies the timepoints (t>>0) to measure GR
X = zeros(length(time)+1,length(dose));
GR = zeros(length(time)+1,length(dose));

[~,dose_id] = find(ismember([O.Dose],dose)==1);

n = 0;

for d = dose_id
    n = n+1;
    for t = 1:length(time)+1
        if t==1
            X(t,:) = O(d).Nlive(1);
        else
            [~,ti] = min(abs(O(d).Tedges-time(t-1))); %index of "O.T" where it is closest to the timepoint specified in "time";
            X(t,n) = O(d).Nlive(ti);
        end
    end
    
end
X_norm = X./X(1,:);
for t = 1:length(time)+1
    GR(t,:) = 2.^(log2(X_norm(t,:))/log2(X_norm(t,1)))-1;
end
end


function [x_CI] = getCI(x,CI)

SEM = nanstd(x,0,1)./sqrt(size(x,1)); % Standard Error
ts = tinv(CI,size(x,1)-1);          % T-Score
x_CI(:,1) = nanmean(x,1) + ts(1).*SEM;      % Lower Confidence Intervals
x_CI(:,2) = nanmean(x,1) + ts(2).*SEM;      % Upper Confidence Intervals
end







