function [DIPmat,DIP_stats] = getDIP_with_stats(SimFolder,simID,time,dose,file_prefix)
% Author:        Natacha Comandante Lou (natacom@umich.edu)
% Copyright 2019, All Rights Reserved
%
% For paper "Phenotype-Based Probabilistic
% Analysis of Heterogeneous Responses to Cancer Drugs and Their Combination
% Efficacy" by N.Comandante-Lou, M.Khaliq, D.Venkat, M.Manikkam,
% M.Fallahi-Sichani
% ------------------------------------------------------------------------
% * Description: This function estimates the drug-induced proliferation (DIP) rates
% (developed by Harris et al. 2016) from simulations and calculates their 
% statistics (mean and confidence interval)
%
% * INPUT:
%   - SimFolder: simulation output folder
%   - simID: simulation index
%   - time: time point to measure DIP (hr)
%   - dose: dose for the condition simulated
%   - file_prefix: prefix for the simulation output file
%   
% * OUTPUT:
% - DIPmat(c).mat(s,d,t): DIP metrics estimated from the s th 
%   simultion of condition c, evaluated at specified dose d and timebin t
% - DIP_stats(c).mat(d,:,t): statistics of DIP rates
%   estimated cross all simulations in DIPmat(c).mat(s,d,t). columns: mean,low CI, high CI

Home  = pwd;
for n = simID
    cd(SimFolder);
    load(strcat(file_prefix,int2str(n),'.mat'));
    cd(Home)
    for j = 1:size(simO,2)
        [~,DIP] = getDIP(simO(:,j),time,dose);
        DIPmat(j).mat(n,:,:) = DIP(2:end,:)';
    end
end
CI = [0.025  0.975]; %confidence interval
for j = 1:length(DIPmat)
    for t = 1:length(time)
        DIP_stats(j).mat(:,1,t) = nanmean(DIPmat(j).mat(:,:,t),1);
        DIP_stats(j).mat(:,2:3,t) = getCI(DIPmat(j).mat(:,:,t),CI);
    end
end



end

function [X,DIP] = getDIP(O,time,dose)
X = zeros(length(time)+1,length(dose));
DIP = zeros(length(time)+1,length(dose));

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
    DIP(t,:) = log(X_norm(t,:))/log(X_norm(t,1));
end


end

function [x_CI] = getCI(x,CI)

SEM = nanstd(x,0,1)./sqrt(size(x,1)); % Standard Error
ts = tinv(CI,size(x,1)-1);          % T-Score
x_CI(:,1) = nanmean(x,1) + ts(1).*SEM;      % Lower Confidence Intervals
x_CI(:,2) = nanmean(x,1) + ts(2).*SEM;      % Upper Confidence Intervals
end

