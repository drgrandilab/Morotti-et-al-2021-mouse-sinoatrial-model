%% Analysis of the properties of a SAN cell: HR and APD vs. Time
function [ap_time cl_array apd_array] = function_SAN_vs_Time_analysis(t,Vm,dVm)

% Detection of MDP
mdp_ind = find( (dVm(1:length(dVm)-1)<0) & (dVm(2:length(dVm))>0) );
mdp_count = length(mdp_ind); % number of events (number of APs is mdp_count-1)

if mdp_count > 3, % at least 3 APs
    ap_time = (t(mdp_ind(1:end-1))+t(mdp_ind(2:end)))/2;
    
    cl_array = zeros(1,length(ap_time));
    apd_array = zeros(1,length(ap_time));
    for iii=1:mdp_count-1,
        ind1 = mdp_ind(iii); ind2 = mdp_ind(iii+1);
        
        CL = t(ind2)-t(ind1); % ms, CL

        Vm_min = Vm(ind1); % mV, MDP
        [Vm_max ind_peak] = max(Vm(ind1:ind2));
        AP_amp = Vm_max-Vm_min; % mV
        f090vm=find(Vm(ind1:ind2)>=Vm_min+10/100*AP_amp);
        APD90 = t(ind1+f090vm(end)-1)-t(ind1+f090vm(1)-1); % ms, APD90
        
        cl_array(iii) = CL;
        apd_array(iii) = APD90;
    end    
else
    ap_time = [0 1];
    cl_array = [0 0];
    apd_array = [0 0];
end