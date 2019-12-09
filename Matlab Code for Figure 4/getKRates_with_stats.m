function [kmat,kmat_stats,log2normNlive_stat,dNmat,tbins] = getKRates_with_stats(Home,SimFolder,simID,dose,dbin,dt,tspan,file_prefix)
% Author:        Natacha Comandante Lou (natacom@umich.edu)
% Copyright 2019, All Rights Reserved
%
% For paper "Phenotype-Based Probabilistic
% Analysis of Heterogeneous Responses to Cancer Drugs and Their Combination
% Efficacy" by N.Comandante-Lou, M.Khaliq, D.Venkat, M.Manikkam,
% M.Fallahi-Sichani
% ------------------------------------------------------------------------
% * Description: This function estimates probabilistic phenotype metrics 
% (k's) from simulations and calculates their statistics 
%
% * INPUT:
%   - Home: home folder path
%   - SimFolder: simulation output folder
%   - simID: simulation index
%   - dose: dose for the conditions simulated
%   - dbin: bin size for k's estimation (hr) 
%   - dt: time resolution of the simulation output (hr)
%   - tspan: start time and end time of the simulation (hr)
%   - file_prefix: prefix for the simulation output file
%   
% * OUTPUT:
%   - kmat(c).death/divi/net(s,t,d): phenotype metrics estimated from the s th 
%     simultion of condition c, evaluated at specified dose d and timebin t
%   - kmat_stats(c).death/divi/net(t,:,d): statistics of phenotype rates
%   estimated cross all simulations in kmat(c).death/divi/net(s,t,d). columns: mean,low CI, high CI
%   - tbins: center time bins corresponding to the estimated phenotype metric values (hr)

bedges = tspan(1):dbin:tspan(end);
tedges = tspan(1):dt:tspan(end);

for n = simID
    cd(SimFolder);
    load(strcat(file_prefix,int2str(n),'.mat'));
    cd(Home)
    for j = 1:size(simO,2)
        [k,dN,tbins] = getKRates(simO(:,j),dose,dbin);
        
        kmat(j).death(n,:,:) = reshape([k(:).death],length(tbins),length(dose));
        kmat(j).divi(n,:,:) = reshape([k(:).division],length(tbins),length(dose));
        kmat(j).net(n,:,:) = reshape([k(:).net],length(tbins),length(dose));
        
        dNmat(j).death(n,:,:) = reshape([dN(:).death],length(tbins),length(dose));
        dNmat(j).divi(n,:,:) = reshape([dN(:).division],length(tbins),length(dose));
        dNmat(j).live(n,:,:) = reshape([dN(:).meanlive],length(tbins),length(dose));
        
        N = reshape([simO(:,j).Nlive],length(tedges),length(dose));
        log2normNlive(j).mat(n,:,:) = log2(N./N(1,:));
    end
end
for j = 1:size(simO,2)
    kmat(j).death(isinf(kmat(j).death))=nan;
    kmat(j).divi(isinf(kmat(j).divi))=nan;
    kmat(j).net(isinf(kmat(j).net))=nan;
    log2normNlive(j).mat(isinf(log2normNlive(j).mat))=nan;
end


CI = [0.025  0.975]; 
for j = 1:length(kmat)
    for d = 1:length(dose)
        kmat_stats(j).death(:,1,d) = nanmean(kmat(j).death(:,:,d),1);
        kmat_stats(j).death(:,2:3,d) = getCI(kmat(j).death(:,:,d),CI);
        kmat_stats(j).divi(:,1,d) = nanmean(kmat(j).divi(:,:,d),1);
        kmat_stats(j).divi(:,2:3,d) = getCI(kmat(j).divi(:,:,d),CI);
        kmat_stats(j).net(:,1,d) = nanmean(kmat(j).net(:,:,d),1);
        kmat_stats(j).net(:,2:3,d) = getCI(kmat(j).net(:,:,d),CI);
        if j < 4
        log2normNlive_stat(j).mat(:,1,d) = nanmean(log2normNlive(j).mat(:,:,d),1);
        log2normNlive_stat(j).mat(:,2:3,d) = getCI(log2normNlive(j).mat(:,:,d),CI);
        end
        
    end
end


end

function [k,dN,tbins] = getKRates(O,dose,dbin)

[~,dose_id] = find(ismember([O.Dose],dose)==1); % convert the selected dose as the dose-index in the simO
for j = 1:length(dose_id)
    dt = mean(diff(O(j).Tedges));
    nt = floor(dbin/dt);
    
    dN(j).death = sum(reshape(O(j).dndeath,nt,[]),1);
    dN(j).division = sum(reshape(O(j).dndivi,nt,[]),1);
    N = reshape(O(j).Nlive(2:end),nt,[]);
    dN(j).meanlive = mean(N,1); % mean number of live cells in each bin
    dN(j).difflive = N(end,:) - N(1,:); % delta number of live cells in each bin
    
    k(j).death = dN(j).death./(dN(j).meanlive*dbin);
    k(j).division = dN(j).division./(dN(j).meanlive*dbin);
    k(j).net = dN(j).difflive./(dN(j).meanlive*dbin);
    
    
end
tbins = O(j).Tedges(floor((nt)/2):nt:end);

end


function [x_CI] = getCI(x,CI) 
%Calculate confidence invervals
SEM = nanstd(x,0,1)./sqrt(size(x,1)); % Standard Error
ts = tinv(CI,size(x,1)-1);          % T-Score
x_CI(:,1) = nanmean(x,1) + ts(1).*SEM;      % Lower Confidence Intervals
x_CI(:,2) = nanmean(x,1) + ts(2).*SEM;      % Upper Confidence Intervals
end

