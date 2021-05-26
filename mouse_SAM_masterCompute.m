%% Mouse SAM model - Master file
% Morotti et al. Intracellular Na+ Modulates Pacemaking Activity in Murine
% Sinoatrial Node Myocytes: An In Silico Analysis. Int. J. Mol. Sci. 2021,
% 22(11), 5645; https://doi.org/10.3390/ijms22115645

clear
close all
clc
%% Parameterization selection
% set model_index = 0 for original parameterization by Kharche et al
% set model_index = 1 for introducing modifications in ICaL, Ito, Isus, and If only
% set model_index = 2 for optimized model

model_index = 2; 
%% Loading initial conditions

if model_index == 0
    load yfin_Kharche % model_index = 0;
    disp('Original Kharche et al. model')
elseif model_index == 1
	load yfin_Kharche_updated_currents % model_index = 1;
    disp('Kharche et al. model w/ updated ICaL, Ito, Isus, and If')
elseif model_index == 2
    load yfin_Kharche_optimized % model_index = 2;
    disp('Optimized model of the murine SAMs')
end

y0n = yfinal;
%% Input parameters
% Flag for Na-clamp
Na_clamp = 0; % [0 for free Na, 1 for Na clamp]
if Na_clamp == 1
    disp('Na clamped')
    %y0n(35) = 12;
end

% Isoproterenol administration
ISO = 0; % (0 for control, 1 for ISO - not used)

% Protocol selection
block_index = 1;
% 0 for no stimulation
% 1 for no stimulation & NKA block at 10 s 
% 2 for no stimulation & NCX block at 10 s 
% 3 for no stimulation & LTCC block at 10 s 
% 4 for no stimulation & NKA/NCX/LTCC modulation at 10 s 

% For protocols 1-3:
block_degree = 0.6; % (0 normal function, 1 full block)

% For protocol 4:
block_array = [0.5 0.5 0]; % differential block for NKA/NCX/LTCC

% Sensitivity analysis parameters
par_SA = ones(1,18); % not used
    
p = [model_index Na_clamp ISO block_index block_degree block_array par_SA];

% Duration
duration = 130e3;
%% Run simulation

disp('Running the simulation...')
options = odeset('RelTol',1e-5,'MaxStep',1);
[t,y] = ode15s(@mouse_SAM_eccODEfile,[0 duration],y0n,options,p);

%% Saving final conditions

% Find new ICs (MDP):
%[~, mdp_index] = min(y(end-3000:end,37));
%yfinal = y(mdp_index+(length(y(:,37))-3001),:);

%save yfin_Kharche yfinal
%save yfin_Kharche_updated_currents yfinal
%save yfin_Kharche_optimized yfinal
%% Plot

figure, set(gcf,'color','w')
subplot(3,1,1),set(gca,'box','off','tickdir','out','fontsize',12)
hold on,plot(t*1e-3,y(:,37)), ylabel('Em (mV)')
subplot(3,1,2),set(gca,'box','off','tickdir','out','fontsize',12)
hold on,plot(t*1e-3,y(:,32)*1e3), ylabel('[Ca]i (uM)')
subplot(3,1,3),set(gca,'box','off','tickdir','out','fontsize',12)
hold on,plot(t*1e-3,y(:,35)), ylabel('[Na]i (mM)')
xlabel('Time (s)')
%% Analysis

currents_flag = 0; % set currents_flag to 1 to calculate and plot ion currents

analysis_flag = 1; % set analysis_flag to 1 to extract (and analyze) AP properties

if currents_flag == 1
    disp('Calculating timecourse for currents and other intermediates...')
    
    tc = t; yc = y; % all elements
   
    currents = calcCurrents(tc,yc,p,'currents');
    % currents = 1/capacitance * [v_dot*capacitance ist ina_ttxs ina_ttxr...
    %         icat ical12 ical13 ih ik1 ikr iks ito isus...
    %         ibna ibca ibk inak inaca Jrel*capacitance];

    dvdt = currents(:,1); 
    i_st = currents(:,2); i_na11 = currents(:,3); i_na15 = currents(:,4);
    i_cat = currents(:,5); i_cal12 = currents(:,6); i_cal13 = currents(:,7);
    i_h = currents(:,8); i_k1 = currents(:,9); i_kr = currents(:,10);
    i_ks = currents(:,11); i_to = currents(:,12); i_sus = currents(:,13);
    i_bna = currents(:,14); i_bca = currents(:,15); i_bk = currents(:,16);
    i_nak = currents(:,17); i_naca = currents(:,18); j_rel = currents(:,19);
    j_up = currents(:,20); i_h_na = currents(:,21); i_b = i_bna + i_bca + i_bk;

    figure, set(gcf,'color','w') 
    subplot(4,2,1),set(gca,'box','off','tickdir','out','fontsize',12)
    hold on,plot(t*1e-3,y(:,37)), ylabel('Em (mV)')
    subplot(4,2,3),set(gca,'box','off','tickdir','out','fontsize',12)
    hold on,plot(t*1e-3,y(:,32)*1e6), ylabel('[Ca]i (nM)')
    subplot(4,2,5),set(gca,'box','off','tickdir','out','fontsize',12)
    hold on,plot(t*1e-3,y(:,34)), ylabel('[Ca]x (mM)')
    plot(t*1e-3,10*y(:,33))
    legend('[Ca]SR-up','10 x [Ca]SR-rel')
    subplot(4,2,7),set(gca,'box','off','tickdir','out','fontsize',12)
    hold on,plot(tc*1e-3,i_h,tc*1e-3,i_naca,tc*1e-3,i_nak)
    ylabel('Other (A/F)')
    legend('If','INaCa','INaK')
    xlabel('Time (s)')
    
    subplot(4,2,2),set(gca,'box','off','tickdir','out','fontsize',12)
    hold on,plot(t*1e-3,y(:,35)), ylabel('[Na]i (mM)')
    subplot(4,2,4),set(gca,'box','off','tickdir','out','fontsize',12)
    hold on,plot(tc*1e-3,i_cal13,tc*1e-3,i_cat,tc*1e-3,i_bca)
    ylabel('Ca (A/F)')
    legend('ICaL','ICaT','IbCa')
    subplot(4,2,6),set(gca,'box','off','tickdir','out','fontsize',12)
    hold on,plot(tc*1e-3,i_na15,tc*1e-3,i_na11,tc*1e-3,i_st,tc*1e-3,i_bna)
    ylabel('Na (A/F)')
    legend('INa1.5','INa1.1','Ist','IbNa')
    subplot(4,2,8),set(gca,'box','off','tickdir','out','fontsize',12)
    hold on,plot(tc*1e-3,i_to,tc*1e-3,i_kr,tc*1e-3,i_k1)
    plot(tc*1e-3,i_sus,tc*1e-3,i_ks,tc*1e-3,i_bk)
    ylabel('K (A/F)')
    legend('Ito','IKr','IK1','Isus','IKs','IbK')
    xlabel('Time (s)')
    
else
    if analysis_flag == 1
        disp('Calculating timecourse for dEm/dt...')
        
        tc = t; yc = y;
        currents = calcCurrents(tc,yc,p,'dVm');
        dvdt = currents(:,1);   
    end
end

if analysis_flag == 1
    time = t; % (ms)
    Vm = y(:,37); % (mV)
    Ca = y(:,32); % (mM)
    Na = y(:,35); % (mM)
    dVm = dvdt; % (mV/ms)

    [ap_time cl_array apd_array] = function_SAN_vs_Time_analysis(t,Vm,dVm);
    
    figure, set(gcf,'color','w')
    subplot(4,1,1),set(gca,'box','off','tickdir','out','fontsize',12)
    hold on,plot(t*1e-3,Vm), ylabel('Em (mV)'), plot_x_lim = xlim;
    subplot(4,1,2),set(gca,'box','off','tickdir','out','fontsize',12)
    hold on,plot(t*1e-3,Ca*1e6), ylabel('[Ca]i (nM)')
    subplot(4,1,3),set(gca,'box','off','tickdir','out','fontsize',12)
    hold on,plot(t*1e-3,Na), ylabel('[Na]i (mM)')
    subplot(4,1,4),set(gca,'box','off','tickdir','out','fontsize',12)
    hold on,plot(ap_time*1e-3,60*1000./cl_array)
    ylabel('FR (bpm)')
    xlim(plot_x_lim)
    xlabel('Time (s)')
    
    Na_HR = zeros(size(ap_time));
    for i = 1:length(ap_time)
        index_array = find(t>ap_time(i));
        Na_HR(i) = y(index_array(1),35);
    end
    figure, set(gcf,'color','w')
    plot(Na_HR,60*1000./cl_array)
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('[Na]i (mM)'), ylabel('FR (bpm)')
    
    % Analysis AP properties - last beat
    disp('Biomarkers analysis:')
    newoutputs = function_SAN_AP_analysis_single_beat(t,Vm,Ca,Na,dVm,0,2);
    % rr_ms rr_Hz rr_bpm (average values over multiple beats)
    % dVm_max -dVm_min Vm_max -Vm_min AP_amp (single beat)
    % -TOP CL APD APD90 APD50 DD eDD eDDR lDD DDR 
    % Ca_min Ca_max Ca_amp Ca_t50 Ca_tau Na_min

    outputs = newoutputs;
    outputs(16) = newoutputs(16)*1e6; % mM -> nM
    outputs(17) = newoutputs(17)*1e6; % mM -> nM
    
    biomarkers = outputs'
    HR_bpm = newoutputs(1)
end   
