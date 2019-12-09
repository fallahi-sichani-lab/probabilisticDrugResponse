function [O] = UnifyTstep(O,tspan,ndt)
%This function converts event times (output from the gillepie algorithm)
%into the number of events occured in each uniformly sized small interval dt
%(dt=tspan/ndt)
% INPUT:
% ndt - number of small interval
% tspan - [simulation start time, simulation end time]
% OUTPUT:
% O.dndeath/dndivi/dnstat - number of death/division/stasis events during each dt
% O.Nlive - number of live cell as a function of time
% O.Tedges - time points (hr) corresponding to the output

[row,col] = size(O);
for j = 1:row
    for k = 1:col
        Tedges = linspace(tspan(1),tspan(2),ndt+1);
        
        dt = mean(diff(Tedges));
        Tedges(1) = Tedges(1) + 1e-5*dt;
        
        [O(j,k).dndivi,~] = histcounts(O(j,k).Count.tdivi(2:end),Tedges);
        [O(j,k).dndeath,~] = histcounts(O(j,k).Count.tdeath(2:end),Tedges);
        [O(j,k).dnstat,~] = histcounts(O(j,k).Count.tstat(2:end),Tedges);
        Xic = O(j,k).X(1,1);
        O(j,k).Nlive = Xic + cumsum(O(j,k).dndivi - O(j,k).dndeath);
        O(j,k).Nlive = [Xic, O(j,k).Nlive];
        O(j,k).Tedges = Tedges;
        
    end
end

