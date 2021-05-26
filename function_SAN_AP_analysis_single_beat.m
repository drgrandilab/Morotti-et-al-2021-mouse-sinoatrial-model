%% Analysis of the properties of a SAN cell (multiple beats)
function newoutputs = function_SAN_AP_analysis_single_beat(t,Vm,Ca,Na,dVm,plot_flag,beat_flag)
% Different analysis based on modifications in experimental data analysis
% by Christian Rickert (November 2016)

% 2018 update: flag "beat_flag" for selection between 1st and last beat
% beat_flag = 1 (analysis of 1st beat), otherwise (analysis of last beat)

% Only the output "rr_bpm" is estimated from the whole simulation,
% independently from the selected beat (beat_flag)

%% Detection of Em Max, mean interval/frequency values
% AP peak
pk_ind = find( (dVm(1:length(dVm)-1)>0) & (dVm(2:length(dVm))<0) );
pk_count = length(pk_ind); % number of events
%pk_freq = 60*1000*pk_count / (t_2(pk_ind(length(pk_ind))) - t_2(pk_ind(1))); % events/min

% Mean peak-to-peak interval and frequency values
rr_bpm = 0;
if pk_count > 3 && (max(Vm)-min(Vm)) > 0.50
    rr_array = t(pk_ind(2:end)) - t(pk_ind(1:end-1));
    rr_ms = mean(rr_array); % ms, peak-to-peak interval (mean value)
    rr_Hz = 1/(rr_ms/1000); % Hz, peak frequency (mean value)
    rr_bpm = 60*rr_Hz; % bpm, peak frequency (mean value)
end

%% Detection of Em Min, SAN properties
% Min diast potential
mdp_ind = find( (dVm(1:length(dVm)-1)<0) & (dVm(2:length(dVm))>0) );
mdp_count = length(mdp_ind); % number of events
%mdp_freq = 60*1000*mdp_count / (t_2(mdp_ind(length(mdp_ind))) - t_2(mdp_ind(1))); % events/min

% Analysis of beat properties
% This function analyzes the action potential (AP) properties based on the 
% membrane potential signal (Vm) and its first derivative (dVm/dt).
% Minimum diastolic potential (MDP) and peak potential (Vmax) were defined as
% the most negative and positive Vm values, respectively.
% Their difference (Vmax-MDP) determines the AP amplitude (APamp). 
% Cycle length was defined as the interval between MDPs of successive APs.
% AP durations at 50 and 90% of repolarization (APD50 and APD90) were calculated
% as the time intervals in which Vm is above MDP+50% APamp and MDP+10% APamp,
% respectively.
% The upstroke velocity (dV_max) was assessed as the maximum value of dVm/dt.
% The repolarization rate (dV_min) was assessed as dVm/dt minimum. We also
% differentiated the minimum dVm/dt values before and after 50% AP repolarization
% (transient repolarization rate, TRR, and maximal repolarization rate, MRR,
% respectively).
% The values of dVm/dt between MDP and the first point used for determining 
% the APD90 interval were fitted (linear regression). The threshold potential
% (THR) was defined as the potential at a time point where dVm/dt exceeds 0.5
% mV/ms above the described linear fit.
% AP duration (APD) was then defined as the interval % between THR and subsequent
% MDP, while distolic duration (DD = CL - APD) as the interval between THR and
% the previous MDP. Diastolic depolarization rate (DDR) was obtained as the
% slope of a linear fit of Vm values between 10 and 50% of the DD interval.
% Late diastolic depolarization rate (lDDR) was obtained as the slope of a
% linear fit of Vm values between 30 and 50% of the DD interval.

if mdp_count > 3 && (max(Vm)-min(Vm)) > 0.50
    if beat_flag == 1   % first beat
        ind1 = mdp_ind(1); ind2 = mdp_ind(2);
    else                % last beat
        ind1 = mdp_ind(end-1); ind2 = mdp_ind(end);
    end
    
    [dVm_max ind_UV] = max(dVm(ind1:ind2)); % mV/ms, upstroke velocity
    [dVm_min ind_RR] = min(dVm(ind1:ind2)); % mV/ms, repolarization rate
    
    %f1=find(t>=t_2(ind1)); t_ind1=f1(1);
    %f2=find(t>=t_2(ind2)); t_ind2=f2(1);

    Vm_min = Vm(ind1); % mV, minimum diastolic potential
    [Vm_max ind_peak] = max(Vm(ind1:ind2));
    AP_amp = Vm_max-Vm_min; % mV
    
    CL = t(ind2)-t(ind1); % ms, cycle length
    
    f090vm=find(Vm(ind1:ind2)>=Vm_min+10/100*AP_amp);
    APD90 = t(ind1+f090vm(end)-1)-t(ind1+f090vm(1)-1); % ms, APD90
    
    f050vm=find(Vm(ind1:ind2)>=Vm_min+50/100*AP_amp);
    APD50 = t(ind1+f050vm(end)-1)-t(ind1+f050vm(1)-1); % ms, APD50
        %TRR = min(dVm(ind1:ind1+f050vm(end)-1)); % mV/ms, TRR
        [late_dVm_min MRR_ind] = min(dVm(ind1+f050vm(end)-1):ind2); % mV/ms, MRR
        
    p = polyfit(t(ind1:ind1+f090vm(1)-1),dVm(ind1:ind1+f090vm(1)-1),1);
    fit_dVm = p(1)*t(ind1:ind2)+p(2);
    f_THR = find(dVm(ind1:ind2)>fit_dVm+0.5);
    if isempty(f_THR) == 0,
        i_THR = ind1+f_THR(1)-1;
        t_THR = t(i_THR);
        THR = Vm(i_THR); % mV, threshold potential

        DD = t_THR - t(ind1); % ms, diastolic duration
        APD = CL - DD;

        EDD = DD/2; % LDD = EDD;
        fDDvm1 = find(t(ind1:ind2)>=t(ind1)+DD*(1/10)); % from 10% DD
        fDDvm2 = find(t(ind1:ind2)<=t(ind1)+DD*(5/10)); % to 50% DD
        q = polyfit(t(ind1+fDDvm1(1)-1:ind1+fDDvm2(end)-1),Vm(ind1+fDDvm1(1)-1:ind1+fDDvm2(end)-1),1);
        fit_Vm = q(1)*t(ind1:ind2)+q(2);
        DDR = q(1); % mV/ms, diastolic depolarization rate

            late_fDDvm1 = find(t(ind1:ind2)>=t(ind1)+DD*(3/10)); % from 30% DD
            late_q = polyfit(t(ind1+late_fDDvm1(1)-1:ind1+fDDvm2(end)-1),Vm(ind1+late_fDDvm1(1)-1:ind1+fDDvm2(end)-1),1);
            fit2_Vm = late_q(1)*t(ind1:ind2)+late_q(2);
            late_DDR = late_q(1); % mV/ms, late diastolic depolarization rate
    else
        THR = 0; % mV, threshold potential
        DD = 0; % ms, diastolic duration
        APD = 0;
        EDD = DD/2; % LDD = EDD;
        DDR = 0; % mV/ms, diastolic depolarization rate
            late_DDR = 0; % mV/ms, late diastolic depolarization rate
    end
    
    % Ca
    Ca_min = min(Ca(ind1:ind2)); % mM, diastolic [Ca]i
    [Ca_max ind_ca_peak] = max(Ca(ind1:ind2)); % mM, systolic [Ca]i
    Ca_amp = Ca_max-Ca_min; % mM, Ca transient amplitude
    
    fCa50 = find(Ca(ind1:ind2)>=Ca_max-0.5*Ca_amp);
    Ca_t50 = t(ind1+fCa50(end)-1)-t(ind1+ind_ca_peak-1); % ms, 50% CaT decay
    
    fCa63 = find(Ca(ind1:ind2)>=Ca_max-0.63*Ca_amp);
    Ca_tau = t(ind1+fCa63(end)-1)-t(ind1+ind_ca_peak-1); % ms, tau CaT decay
    
    % Na
    %Na_conc = Na(ind1); % mM, [Na]i
    Na_min = min(Na(ind1:ind2)); % mM, [Na]i
    
	% Save AP signal
    %AP_ctrl_y = Vm(ind1:ind2);
    %AP_ctrl_x = t(ind1:ind2)-t(ind1);
    %figure(105),plot(AP_ctrl_x,AP_ctrl_y)
    %save AP_opt_model AP_ctrl_x AP_ctrl_y
    
    % Plot figure (based on input plot_flag)
    if plot_flag == 1 % Figure
        
        figure(101),set(gcf,'color','w')
        subplot(2,1,1),hold on,set(gca,'box','off','tickdir','out','fontsize',12)
        plot(t(ind1:ind2)-t(ind1),Vm(ind1:ind2)) % t0 @ t_in
        ylabel('Em (mV)')
        subplot(2,1,2),hold on,set(gca,'box','off','tickdir','out','fontsize',12)
        plot(t(ind1:ind2)-t(ind1),dVm(ind1:ind2)) % t0 @ t_in
        ylabel('dEm/dt (mV/ms)'), xlabel('Time (ms)')

        subplot(2,1,1)
        plot(t(ind1+ind_peak-1)-t(ind1),Vm_max,'or')
        plot([0 CL],[Vm_min Vm_min],'or')
        plot([0 CL],(Vm_max+25)*[1 1],'k')
        plot([CL+10 CL+10],[Vm_min Vm_max],'k')
        plot([0 CL],Vm_min+0.50*AP_amp*[1 1],'k:',[0 CL],Vm_min+0.10*AP_amp*[1 1],'k:')
        plot([t(ind1+f050vm(1)-1) t(ind1+f050vm(end)-1)]-t(ind1),Vm_min+0.50*AP_amp*[1 1],'k')
        plot([t(ind1+f090vm(1)-1) t(ind1+f090vm(end)-1)]-t(ind1),Vm_min+0.10*AP_amp*[1 1],'k')
        plot((t(ind1+f050vm(end)-1)-t(ind1))*[1 1],[Vm_min Vm_max],'k:')
        plot((t(ind1+f090vm(1)-1)-t(ind1))*[1 1],[Vm_min Vm_max],'k:')
        if isempty(f_THR) == 0,
            plot(t_THR-t(ind1),THR,'om')
            plot((t_THR-t(ind1))*[1 1],[Vm_min Vm_max],'m:')
            plot([t_THR-t(ind1) t_THR-t(ind1)+APD],(Vm_max+15)*[1 1],'k')
            plot([0 EDD],(Vm_max+15)*[1 1],'k')
            plot_DDR = 0;
            if plot_DDR == 1,
                plot(t(ind1+fDDvm1(1)-1:ind1+fDDvm2(end)-1)-t(ind1),Vm(ind1+fDDvm1(1)-1:ind1+fDDvm2(end)-1),'c')
                plot((t(ind1+fDDvm1(1)-1)-t(ind1))*[1 1],[Vm_min Vm_max],'c--')
                plot((t(ind1+fDDvm2(end)-1)-t(ind1))*[1 1],[Vm_min Vm_max],'c--')
                plot(t(ind1:ind2)-t(ind1),fit_Vm,'c--');
            end
            plot_lDDR = 1;
            if plot_lDDR == 1,
                plot(t(ind1+late_fDDvm1(1)-1:ind1+fDDvm2(end)-1)-t(ind1),Vm(ind1+late_fDDvm1(1)-1:ind1+fDDvm2(end)-1),'g')
                plot((t(ind1+late_fDDvm1(1)-1)-t(ind1))*[1 1],[Vm_min Vm_max],'g:')
                plot((t(ind1+fDDvm2(end)-1)-t(ind1))*[1 1],[Vm_min Vm_max],'g:')
                plot(t(ind1:ind2)-t(ind1),fit2_Vm,'g--');
            end
        end
        xlim([-10 CL+15])

        subplot(2,1,2)
        plot([t(ind1+ind_RR-1) t(ind1+ind_UV-1)]-t(ind1),[dVm_min dVm_max],'or')
        plot([CL+10 CL+10],[dVm_min dVm_max],'k')
        plot((t(ind1+f090vm(1)-1)-t(ind1))*[1 1],[dVm_min dVm_max],'k:')
        plot((t(ind1+f050vm(end)-1)-t(ind1))*[1 1],[dVm_min dVm_max],'k:')
        plot(t(ind1+f050vm(end)-1+(MRR_ind-1))-t(ind1),late_dVm_min,'*k')
        plot(t(ind1:ind1+f090vm(1)-1)-t(ind1),dVm(ind1:ind1+f090vm(1)-1),'m')
        plot(t(ind1:ind2)-t(ind1),fit_dVm,'m--');
        if isempty(f_THR) == 0,
            plot((t_THR-t(ind1))*[1 1],[dVm_min dVm_max],'m:')
        end
        xlim([-10 CL+15])
        plot(xlim,[0 0],'k:')

        if isempty(f_THR) == 0
            figure(102),set(gcf,'color','w')
            subplot(2,1,1),hold on,set(gca,'box','off','tickdir','out','fontsize',12)
            plot(t(ind1:ind2)-t(i_THR),Vm(ind1:ind2)) % t0 @ THR
            ylabel('Em (mV)')
            subplot(2,1,2),hold on,set(gca,'box','off','tickdir','out','fontsize',12)
            plot(t(ind1:ind2)-t(i_THR),dVm(ind1:ind2)) % t0 @ THR
            ylabel('dEm/dt (mV/ms)'), xlabel('Time (ms)')
        end
        
    end
    
end

%% Collection of all the outputs
% rr_ms rr_Hz rr_bpm (average values over multiple beats)
% dVm_max -dVm_min Vm_max -Vm_min AP_amp (single beat)
% -TOP CL APD APD90 APD50 DD eDD eDDR lDD DDR 
% Ca_min Ca_max Ca_amp Ca_t50 Ca_tau Na_min

if mdp_count > 3 && pk_count > 3 && (max(Vm)-min(Vm)) > 0.50
    newoutputs = [rr_bpm dVm_max -dVm_min -Vm_min AP_amp...
        -THR APD APD90 APD50 CL...
        DD EDD DDR late_DDR -late_dVm_min...
        Ca_min Ca_amp Ca_t50 Ca_tau Na_min];
else
    newoutputs = zeros(1,20);
        newoutputs(4) = -min(Vm);
        newoutputs(16) = min(Ca);
        newoutputs(20) = min(Na);
end