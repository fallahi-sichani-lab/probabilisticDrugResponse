function [O] = UnifyTstep_Combo(O,tspan,ndt)
%This function converts event times (output from the gillepie algorithm)
%into the number of events occured in each uniformly sized small interval dt
%(dt=tspan/ndt)
% INPUT:
% ndt - number of small interval
% tspan - [simulation start time, simulation end time]
% OUTPUT:
% O.dndeath/dndivi - number of death or division events during each dt
% O.Nlive - number of live cell as a function of time
% O.Tedges - time points (hr) corresponding to the output
% O.dndeath_sub - number death events during each dt, induced by drug A
% (1st column), and drug B (2nd column)
[ndose,ncase] = size(O);
Tedges = linspace(tspan(1),tspan(2),ndt+1);
dt = mean(diff(Tedges));
Tedges(1) = Tedges(1) + 1e-5*dt;
for j = 1:ncase
    for d = 1:ndose
        [O(d,j).dndivi,~] = histcounts(O(d,j).Count.tdivi(2:end),Tedges);
        if isstruct(O(d,j).Count.tdeath)
            for k = 1:length(O(d,j).Count.tdeath)
                [O(d,j).dndeath_sub(k,:),~] = histcounts(O(d,j).Count.tdeath(k).deaths(2:end),Tedges);
                
                %O(d,j).dndeath_sub are deaths from individual drugs
            end
            O(d,j).dndeath = sum(O(d,j).dndeath_sub,1);
        else
            [O(d,j).dndeath,~] = histcounts(O(d,j).Count.tdeath(2:end),Tedges);
        end
        
        Xic = O(d,j).X(1,1);
        O(d,j).Nlive = Xic + cumsum(O(d,j).dndivi - O(d,j).dndeath);
        O(d,j).Nlive = [Xic, O(d,j).Nlive];
        O(d,j).Tedges = Tedges;
    end
    
end

