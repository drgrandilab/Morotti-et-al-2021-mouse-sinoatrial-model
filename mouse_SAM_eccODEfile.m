%% Mouse SAM model - ODE file
% Morotti et al. Intracellular Na+ Modulates Pacemaking Activity in Murine
% Sinoatrial Node Myocytes: An In Silico Analysis. Int. J. Mol. Sci. 2021,
% 22(11), 5645; https://doi.org/10.3390/ijms22115645

function output = mouse_SAM_eccODEfile(t,y,p,runType)

ydot = zeros(size(y));
%% Input parameters
% Model index:
% 0 for Kharche model,
% 1 for updated currents,
% 2 for updated currents, & Optimization

m_ind = p(1); 

Na_clamp = p(2); % flag Na-clamp (0 for free Na, 1 for Na clamp)

ISO = p(3); % flag ISO (not used)

block_index = p(4);
% 1 for no stimulation & NKA block at 10 s 
% 2 for no stimulation & NCX block at 10 s 
% 3 for no stimulation & LTCC block at 10 s 
% 4 for no stimulation & NKA/NCX/LTCC modulation at 10 s 

block_input = p(5); % for protocol parameter (V-clamp/NKA block)

block_array = p(6:8); % differential block for NKA/NCX/LTCC

par_SA = p(9:end); % for sensitivity analysis
% 1) gst 2) gna_ttxs 3) gna_ttxr 4) gcat 5) gcal12 6) gcal13 
% 7) gh 8) gk1 9) gkr 10) gks 11) gto 12) gsus
% 13) gbna 14) gbca 15) inakmax 16) kNaCa 17) ks 18) Pup

%% NKA/NCX/LTCC

inakmax_multiplier = 1; 
inaca_multiplier = 1; 
ical_multiplier = 1;

if block_index == 1
    NKA_block = block_input; % 0.50 or 0.75%
    inakmax_multiplier = (1-(t>9.95e3)*NKA_block); % NKA block after 10 s
    %inakmax_multiplier = (1-(t>9.95e3)*(t<39.95e3)*NKA_block); % transient NKA block
end

if block_index == 2
    NCX_block = block_input; % 0.50 or 0.75%
    inaca_multiplier = (1-(t>9.95e3)*NCX_block); % NCX block after 10 s
    %inaca_multiplier = (1-(t>9.95e3)*(t<40e3)*NCX_block); % transient NCX block
end

if block_index == 3
    LTCC_block = block_input; % 0.50 or 0.75%
    ical_multiplier = (1-(t>9.95e3)*LTCC_block); % LTCC block after 10 s
    %ical_multiplier = (1-(t>9.95e3)*(t<39.95e3)*LTCC_block); % transient LTCC block
end

if block_index == 4
    NKA_block = block_array(1); % 0.50 or 0.75%
    inakmax_multiplier = (1-(t>9.95e3)*NKA_block); % NKA block after 10 s
    NCX_block = block_array(2); % 0.50 or 0.75%
    inaca_multiplier = (1-(t>9.95e3)*NCX_block); % NCX block after 10 s
    LTCC_block = block_array(3); % 0.50 or 0.75%
    ical_multiplier = (1-(t>9.95e3)*LTCC_block); % LTCC block after 10 s
end
%% State variables

dst = y(1);
fst = y(2);
dt = y(3);
ft = y(4);
ikr_act = y(5);
ikr_inact = y(6);
iks_act = y(7);
dl13 = y(8);
fl13 = y(9);
dl12 = y(10);
fl12 = y(11);
fca = y(12);
m_ttxs = y(13);
h_ttxs = y(14);
j_ttxs = y(15);
m_ttxr = y(16);
h_ttxr = y(17);
j_ttxr = y(18);
y_1_2 = y(19);
q = y(20);
r = y(21);
resting = y(22);
open = y(23);
inactivated = y(24);
Ftc = y(25);
Ftmc = y(26);
Ftmm = y(27);
Fcms = y(28);
Fcmi = y(29);
Fcq = y(30);
casub = y(31);
cai = y(32);
carel = y(33);
caup = y(34);
nai = y(35);
ki = y(36);
v = y(37);
yb_1_2 = y(38);
%% Model Parameters - Kharche et al model

R = 8.314472; % [J/mol*K]
T = 310.5; % [K] 
F = 96.4845; % [C/kmol]
FRT = F/(R*T);

capacitance = 0.025; % [nF]
vsub = 0.03328117;
vi = 1.34671883;
vrel = 0.0036;
vup = 0.0348;

Mgi = 2.5; % [mM]
nao = 140; % [mM]
cao = 1.8; % [mM]
ko = 5.4; % [mM]

gst=0.00006;
eist=17.0;
gbna=0.0001215;
gbca=0.000015;
gbk =0.0000025;
gk1=0.229*0.0039228*0.9;
gks=0.000299;
ecal=47.0;
kmfca=0.00035;
alpha_fca=0.021;
ecat=45.0;
enattxr=41.5761;
gsus=0.00039060;
inakmax=inakmax_multiplier*1.85*0.077;
kmnap=14.0;
kmkp=1.4;
K1ni=395.3;
K1no=1628;
K2ni=2.289;
K2no=561.4;
K3ni=26.44;
K3no=4.663;
Kci=0.0207;
Kco=3.663;
Kcni=26.44;
Qci=0.1369;
Qco=0.0;
Qn=0.4315;
tdifca=0.04;
Ttr=40.0;

% Buffer
ConcTC=0.031;
kfTC=88.8;
kbTC=0.446;
ConcTMC=0.062;
kfTMC=237.7;
kbTMC=0.00751;
kfTMM=2.277;
kbTMM=0.751;
ConcCM=0.045;
kfCM=237.7;
kbCM=0.542;
ConcCQ=10.0;
kfCQ=0.534;
kbCQ=0.445;

koca=10.0;
kom=0.06;
kica=0.5;
kim=0.005;
eca50sr=0.45;
maxsr=15.0;
minsr=1.0;
hsrr=2.5;
pumphill=2.0;

gna_ttxs = 0.1*5.925e-05; % INa1.1
gna_ttxr = 0.1*5.925e-05; % INa1.5
v_cal = 0;
gcal12 = 0.0010*4.0*1.5;
gcal13 = 0.0030*4.0*1.5;
gcat = 0.75*0.01862; 
gh = 0.0057; 
slope_gh = 16.3;
vhalf_gh = 106.8;
gkr = 1*(0.8*0.002955);
gto = 0.00492;
kNaCa = inaca_multiplier*5.5;
Pup = 0.04;
ks = 1300000;
pumpkmf = 0.00008; % c code
pumpkmr = 4.5; % c code

% flag If model
new_If_flag = 0;
new_If_ISO_flag = 0;

if ISO == 1 % only affects If (see Peteres et al PNAS 2021)
    new_If_ISO_flag = 1;
end
%% Model Parameters - Modified
% m_ind = 0 -> Kharche

if m_ind == 1 % -> Updated currents
    gto = gto*2.5;
    gsus = gsus*2.5;
    
    gcal12 = gcal12*0;
    gcal13 = gcal13*2/3;
    v_cal = 7;
    
    new_If_flag = 1;
    slope_gh = 11.66; vhalf_gh = 103.2; % low (nM) ISO
    gh = gh*3;
end

if m_ind == 2 % -> Updated currents & Optimization
    gto = gto*2.5;
    gsus = gsus*2.5;
    
    gcal12 = gcal12*0;
    gcal13 = gcal13*2/3;
    v_cal = 7;
    
    new_If_flag = 1;
    slope_gh = 11.66; vhalf_gh = 103.2; % low (nM) ISO
    gh = gh*3;
    
    % Optimization
    opt_factors = [0.993562274176361,1.07970775119823,1.32358104543159,2.10193703572732,1,1.61395022486174,1.45398234403856,0.901694268755210,2.23421113694337,0.907645260031040,1.23102830269302,0.861664782770364,1.34506700840916,1.06876306706288,1.53761453695508,0.822428341648628,1.29297095679933,1.08451353877899];
    par_SA = par_SA.*opt_factors;
end
%% Model equations

ena = (R*T/F)*log(nao/nai);
ek  = (R*T/F)*log(ko/ki);
eks = (R*T/F)*log((ko+0.12*nao)/(ki+0.12*nai));
eca = (R*T/(2*F))*log(cao/casub);
%% Ist - Sustained inward Na current **************************************

qa = 1.0/(1.0 + exp(-(v+67.0)/5.0));
alphaqa = 1.0/(0.15*exp(-(v)/11.0)+0.2*exp(-(v)/700.0));
betaqa = 1.0/(16.0*exp((v)/8.0)+15.0*exp((v)/50.0));
tauqa = 1.0/(alphaqa + betaqa);
alphaqi = 0.15*1.0/(3100.0*exp((v+10.0)/13.0)+700.3*exp((v+10.0)/70.0));
betaqi =  0.15*1.0/(95.7*exp(-(v+10.0)/10.0) + 50.0*exp(-(v+10.0)/700.0)) + 0.000229/(1+exp(-(v+10.0)/5.0));
qi = alphaqi/(alphaqi + betaqi);
tauqi = 1.0/(alphaqi + betaqi);
dst_dot = (qa-dst)/tauqa;
fst_dot = (qi-fst)/tauqi;
ist = par_SA(1)*gst*dst*fst*(v - eist);
%% INa - Na channel isoforms Nav1.1/Nav1.5 currents ***********************

fna = (9.52e-02*exp(-6.3e-2*(v+34.4))/(1+1.66*exp(-0.225*(v+63.7))))+8.69e-2; 
m3_inf_ttxr = 1.0/(1.0 + exp(-(v+45.213705)/7.219547));
h_inf_ttxr = 1.0/(1.0 + exp(-(v+62.578120 )/(-6.084036)));
m3_inf_ttxs = 1.0/(1.0 + exp(-(v+36.097331-5.0)/5.0));
h_inf_ttxs = 1.0/(1.0 + exp((v+56.0)/3.0));
m_inf_ttxr = m3_inf_ttxr^0.333;
m_inf_ttxs = m3_inf_ttxs^0.333;
tau_m = 1000.0*((0.6247e-03/(0.832*exp(-0.335*(v+56.7))+0.627*exp(0.082*(v+65.01))))+0.0000492);
tau_h = 1000.0*(((3.717e-06*exp(-0.2815*(v+17.11)))/(1+0.003732*exp(-0.3426*(v + 37.76))))+0.0005977);
tau_j = 1000.0*(((0.00000003186*exp(-0.6219*(v+18.8)))/(1+0.00007189*exp(-0.6683*(v+34.07))))+0.003556);
m_ttxs_dot = (m_inf_ttxs - m_ttxs)/tau_m;
h_ttxs_dot = (h_inf_ttxs - h_ttxs)/tau_h;
j_ttxs_dot = (h_inf_ttxs - j_ttxs)/tau_j;
hs = (1.0-fna)*h_ttxs+fna*j_ttxs;
tau_mr = 1000.0*((0.6247e-03/(0.832*exp(-0.335*(v+56.7))+0.627*exp(0.082*(v+65.01))))+0.0000492);
tau_hr = 1000.0*(((3.717e-06*exp(-0.2815*(v+17.11)))/(1+0.003732*exp(-0.3426*(v + 37.76))))+0.0005977);
tau_jr = 1000.0*(((0.00000003186*exp(-0.6219*(v+18.8)))/(1+0.00007189*exp(-0.6683*(v+34.07))))+0.003556);
m_ttxr_dot = (m_inf_ttxr - m_ttxr)/tau_mr;
h_ttxr_dot = (h_inf_ttxr - h_ttxr)/tau_hr;
j_ttxr_dot = (h_inf_ttxr - j_ttxr)/tau_jr;
hsr = (1.0-fna)*h_ttxr+fna*j_ttxr;
if abs(v)>0.005
    f_ina_ttxs= m_ttxs*m_ttxs*m_ttxs*hs*nao*(F*F/(R*T))*((exp((v-ena)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
    f_ina_ttxr = m_ttxr*m_ttxr*m_ttxr*hsr*nao*(F*F/(R*T))*((exp((v-enattxr)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
    f_ina_ttxs= m_ttxs*m_ttxs*m_ttxs*hs*nao*F*((exp((v-ena)*F/(R*T))-1.0));
    f_ina_ttxr = m_ttxr*m_ttxr*m_ttxr*hsr*nao*F*((exp((v-enattxr)*F/(R*T))-1.0));
end
ina_ttxs= par_SA(2)*gna_ttxs*f_ina_ttxs;
ina_ttxr = par_SA(3)*gna_ttxr*f_ina_ttxr;
%% ICaT Cav3.1 - T-type Ca channel isoform Cav3.1 current *****************

v_cat = 0;
tau_dt = 1.0/(1.068*exp((v + 26.3 + v_cat)/30.0) + 1.068*exp(-(v + 26.3 + v_cat)/30.0));
dt_inf = 1.0/(1.0+exp(-(v + 26.0 + v_cat)/6.0));
dt_dot = (dt_inf - dt)/tau_dt;
tau_ft = 1.0/(0.0153*exp(-(v+61.7 + v_cat)/83.3)+0.015*exp((v+61.7 + v_cat)/15.38));
ft_inf = 1.0/(1.0+exp((v + 61.7 + v_cat)/5.6));
ft_dot = (ft_inf - ft)/tau_ft;
icat = par_SA(4)*gcat*ft*dt*(v - ecat); 
%% ICaL - L-type Ca channel isoforms Cav1.2/Cav1.3 currents ***************

v_cal_shift = 0;
if abs(v+v_cal+v_cal_shift)<=0.001
    alpha_dl = -28.39*(v+35.0+v_cal+v_cal_shift)/(exp(-(v+35.0+v_cal+v_cal_shift)/2.5)-1.0)+408.173;
elseif abs(v+35.0+v_cal+v_cal_shift)<=0.001
    alpha_dl = 70.975-84.9*(v+v_cal+v_cal_shift)/(exp(-0.208*(v+v_cal+v_cal_shift))-1.0);
else
    alpha_dl = -28.39*(v+35.0+v_cal+v_cal_shift)/(exp(-(v+35.0+v_cal+v_cal_shift)/2.5)-1.0)-84.9*(v+v_cal+v_cal_shift)/(exp(-0.208*(v+v_cal+v_cal_shift))-1.0);
end
if abs(v-5.0+v_cal+v_cal_shift)<=0.001
    beta_dl = 28.575;
else
    beta_dl = 11.43*(v-5.0+v_cal+v_cal_shift)/(exp(0.4*(v-5.0+v_cal+v_cal_shift))-1.0);
end
tau_dl = 2000.0/(alpha_dl+beta_dl);
tau_fl = (7.4 + 45.77*exp(-0.5*(v+28.1+v_cal+v_cal_shift)*(v+28.1+v_cal+v_cal_shift)/(11*11)));
dl13_inf = 1.0/(1+exp(-(v+13.5+v_cal+v_cal_shift)/6.0));
fl13_inf = 1.0/(1+exp((v+35.0+v_cal+v_cal_shift)/7.3));
dl13_dot = (dl13_inf - dl13)/tau_dl;
fl13_dot = (fl13_inf - fl13)/tau_fl;
dl12_inf = 1.0/(1+exp(-(v+3.0+v_cal+v_cal_shift)/5.0));
fl12_inf = 1.0/(1+exp((v+36.0+v_cal+v_cal_shift)/4.6));
dl12_dot = (dl12_inf - dl12)/tau_dl;
fl12_dot = (fl12_inf - fl12)/tau_fl;

casub_fca = casub;
fca_inf = kmfca/(kmfca+casub_fca);
taufca = fca_inf/alpha_fca;
fca_dot = (fca_inf - fca)/taufca;

ical12 = ical_multiplier*par_SA(5)*gcal12*fl12*dl12*fca*(v-ecal); % set to 0 in the optimized version
ical13 = ical_multiplier*par_SA(6)*gcal13*fl13*dl13*fca*(v-ecal);
%% If - Hyperpolarization-activated current *******************************

if new_If_flag == 1
    ih_pK = 0.66; % ih_PNa = 0.34;
    % calculated as: ih_pK = (Erev_exp-ena)/(ek-ena), with Erev_exp = -35.4;
    % from: 0 = ih_pK*(Erev_exp-ek) + (1-ih_pK)*(Erev_exp-ena)
    
    % low ISO
    % time constant - slow component
    tau_y_1_2 = 70/(exp(-(v+300)*0.015)+ exp((v-200)*0.015));
    % time constant - fast component
    tau_yb_1_2 = 0.7/(exp(-(v+340)*0.025)+ exp((v-160)*0.025));
    % ratio slow/fast
    ratio_SonF = 0.015/(exp(-(v+300)*0.025)+exp((v-200)*0.025));
        
    if new_If_ISO_flag == 1 % High ISO
        tau_y_1_2 = 10/(exp(-(v+285)*0.025)+ exp((v-215)*0.025)); % high (uM) ISO
        tau_yb_1_2 = 0.035/(exp(-(v+330)*0.040)+ exp((v-170)*0.040)); % high (uM) ISO
        ratio_SonF = 0.003/(exp(-(v+330)*0.030)+ exp((v-170)*0.030)); % high (uM) ISO
        slope_gh = 9.15; vhalf_gh = 88.5; % high (uM) ISO
    end

    % activation
    y_inf = 1.0/(1.0 + exp((v+vhalf_gh)/slope_gh));
        
    y_1_2_dot = (y_inf - y_1_2)/tau_y_1_2; % slow
    yb_1_2_dot = (y_inf - yb_1_2)/tau_yb_1_2; % fast
    
    ih_act = ((ratio_SonF/(ratio_SonF+1)*y_1_2)+(1/(ratio_SonF+1)*yb_1_2));
    
    ihk  = ih_pK*par_SA(7)*gh*ih_act*(v - ek);
    ihna = (1-ih_pK)*par_SA(7)*gh*ih_act*(v - ena);
    
    ih = (ihk + ihna);
else    
    ih_pK = 0.6167; % ih_PNa = 0.3833;
    y_inf = 1.0/(1.0 + exp((v+vhalf_gh)/slope_gh));
    tau_y_1_2 = 1.5049/(exp(-(v+590.3)*0.01094)+ exp((v-85.1)/17.2));

    y_1_2_dot = (y_inf - y_1_2)/tau_y_1_2;
    yb_1_2_dot = 0;

    ihk  = ih_pK*par_SA(7)*gh*y_1_2*(v - ek);
    ihna = (1-ih_pK)*par_SA(7)*gh*y_1_2*(v - ena);
    
    ih = (ihk + ihna);
end
%% IK1 - Time-independent K current ***************************************

xk1inf = 1.0/(1.0 + exp(0.070727*(v - ek)));
ik1 = par_SA(8)*gk1*xk1inf*(ko/(ko + 0.228880))*(v - ek);
%% IKr - Rapid delayed rectifying K current *******************************

v_kr = 0;
ikr_act_inf = 1.0/(1.0 + exp(-(v+21.173694+v_kr)/9.757086 *1));
tau_ikr_act = 1*    ( 0.699821/(0.003596*exp((v)/15.339290) + 0.000177*exp(-(v)/25.868423)) );
ikr_act_dot = (ikr_act_inf-ikr_act)/tau_ikr_act;     
ikr_inact_inf = 1.0/(1.0 + exp((v+20.758474-4.0    +40*0)/(19.0)));
tau_ikr_inact = 1*(0.2+0.9*1.0/(0.1*exp(v/54.645)+0.656*exp(v/106.157)));
ikr_inact_dot = (ikr_inact_inf - ikr_inact)/tau_ikr_inact;
ikr = par_SA(9)*gkr*ikr_act*ikr_inact*(v - ek);
%% IKs - Slow delayed rectifying K current ********************************

iks_act_inf = 1.0/(1.0 + exp(-(v-20.876040)/11.852723));
tau_iks_act = 1000.0/(13.097938/(1.0 + exp(-(v-48.910584)/10.630272)) + exp(-(v)/35.316539));
iks_act_dot = (iks_act_inf - iks_act)/tau_iks_act;
iks = par_SA(10)*gks*iks_act*iks_act*(v - eks);
%% Ito - Transient component of 4-AP-sensitive currents *******************

q_inf = 1.0/(1.0+exp((v+49.0)/13.0));
tau_q = (6.06 + 39.102/(0.57*exp(-0.08*(v+44.0))+0.065*exp(0.1*(v+45.93))))/0.67; 
q_dot = (q_inf-q)/tau_q;
r_inf = 1.0/(1.0+exp(-(v-19.3)/15.0));
tau_r = (2.75+14.40516/(1.037*exp(0.09*(v+30.61))+0.369*exp(-0.12*(v+23.84))))/0.303;
r_dot = (r_inf-r)/tau_r;
ito = par_SA(11)*gto*q*r*(v-ek);
%% Isus - Sustained component of 4-AP-sensitive currents ******************

isus = par_SA(12)*gsus*r*(v-ek);
%% Ib - Background Na, Ca and K currents **********************************

ibna = par_SA(13)*gbna*(v - ena);
ibca = par_SA(14)*gbca*(v - eca);
ibk = gbk*(v - ek);
ib = (ibna + ibca + ibk);
%% INaK - Na-K pump current ***********************************************

inak = par_SA(15)*inakmax*((ko^1.2)/((kmkp^1.2)+(ko^1.2)))*((nai^1.3)/((kmnap^1.3)+(nai^1.3)))/(1.0+exp(-(v-ena+120.0)/30.0));
%% INaCa - Na/Ca exchanger current ****************************************

di=1+(casub/Kci)*(1+exp(-Qci*v*FRT)+nai/Kcni)+(nai/K1ni)*(1+(nai/K2ni)*(1+nai/K3ni));
doo=1+(cao/Kco)*(1+exp(Qco*v*FRT))+(nao/K1no)*(1+(nao/K2no)*(1+nao/K3no));
k43=nai/(K3ni+nai);
k12=(casub/Kci)*exp(-Qci*v*FRT)/di;
k14=(nai/K1ni)*(nai/K2ni)*(1+nai/K3ni)*exp(Qn*v*FRT/2.0)/di;
k41=exp(-Qn*v*FRT/2.0);
k34=nao/(K3no+nao);
k21=(cao/Kco)*exp(Qco*v*FRT)/doo;
k23=(nao/K1no)*(nao/K2no)*(1+nao/K3no)*exp(-Qn*v*FRT/2.0)/doo;
k32=exp(Qn*v*FRT/2);
x1=k34*k41*(k23+k21)+k21*k32*(k43+k41);
x2=k43*k32*(k14+k12)+k41*k12*(k34+k32);
x3=k43*k14*(k23+k21)+k12*k23*(k43+k41);
x4=k34*k23*(k14+k12)+k21*k14*(k34+k32);

inaca = par_SA(16)*kNaCa*(k21*x2-k12*x1)/(x1+x2+x3+x4);
%% SR Ca fluxes ***********************************************************

kcasr = maxsr - (maxsr - minsr)/(1.0 + (eca50sr/carel)^hsrr);
kosrca = koca/kcasr;
kisrca = kica*kcasr;
resting_inactivated = 1-resting-open-inactivated;
resting_dot = (kim*resting_inactivated - kisrca*casub*resting - kosrca*casub*casub*resting + kom*open);
open_dot = (kosrca*casub*casub*resting - kom*open - kisrca*casub*open + kim*inactivated);
inactivated_dot = (kisrca*casub*open - kim*inactivated - kom*inactivated + kosrca*casub*casub*resting_inactivated);

Jrel = par_SA(17)*ks*open*(carel - casub);
Jup = par_SA(18)*Pup*((cai/pumpkmf)^pumphill - (caup/pumpkmr)^pumphill)/(1.0 + (cai/pumpkmf)^pumphill + (caup/pumpkmr)^pumphill);
Jtr  = (caup - carel)/Ttr;
%% Ca Buffering ***********************************************************

Ftc_dot = (kfTC*cai*(1.0-Ftc)-kbTC*Ftc);
Ftmc_dot = (kfTMC*cai*(1.0-Ftmc-Ftmm)-kbTMC*Ftmc);
Ftmm_dot = (kfTMM*Mgi*(1.0-Ftmc-Ftmm)-kbTMM*Ftmm);
Fcms_dot = (kfCM*casub*(1.0-Fcms)-kbCM*Fcms);
Fcmi_dot = (kfCM*cai*(1.0-Fcmi)-kbCM*Fcmi);
Fcq_dot = (kfCQ*carel*(1.0-Fcq)-kbCQ*Fcq);
%% Ca, Na and K concentrations ********************************************

ca_flux = (ical12+ical13+icat-2.0*inaca+ibca)/(2.0*F);
Jcadif = (casub - cai)/tdifca;

casub_dot = ((-ca_flux+Jrel*vrel)/vsub-Jcadif - ConcCM*Fcms_dot); % ConcCM
cai_dot = ((Jcadif*vsub-Jup*vup)/vi - (ConcCM*Fcmi_dot + ConcTC*Ftc_dot + ConcTMC*Ftmc_dot)); 
carel_dot = (Jtr - Jrel - ConcCQ*Fcq_dot);
caup_dot = (Jup-Jtr*vrel/vup);

nai_tot = ihna+ina_ttxr+ina_ttxs+3.0*inak+3.0*inaca+ist+ibna;
nai_dot = (-nai_tot)/(F*vi);

ki_tot = ihk+iks+ikr+ik1+ibk-2.0*inak+isus+ito;
ki_dot = (-ki_tot)/(F*vi);
%% Membrane potential *****************************************************

I_app = 0;

total_current = ih+ina_ttxr+ina_ttxs+ical12+ical13+iks+ikr+ik1+ist+ib+icat+inak+isus+inaca+ito;
v_dot = -(total_current-I_app)/capacitance;
%% Outputs

ydot(1) = dst_dot;
ydot(2) = fst_dot;
ydot(3) = dt_dot;
ydot(4) = ft_dot;
ydot(5) = ikr_act_dot;
ydot(6) = ikr_inact_dot;
ydot(7) = iks_act_dot;
ydot(8) = dl13_dot;
ydot(9) = fl13_dot;
ydot(10) = dl12_dot;
ydot(11) = fl12_dot;
ydot(12) = fca_dot;
ydot(13) = m_ttxs_dot;
ydot(14) = h_ttxs_dot;
ydot(15) = j_ttxs_dot;
ydot(16) = m_ttxr_dot;
ydot(17) = h_ttxr_dot;
ydot(18) = j_ttxr_dot;
ydot(19) = y_1_2_dot;
ydot(20) = q_dot;
ydot(21) = r_dot;
ydot(22) = resting_dot;
ydot(23) = open_dot;
ydot(24) = inactivated_dot;
ydot(25) = Ftc_dot;
ydot(26) = Ftmc_dot;
ydot(27) = Ftmm_dot;
ydot(28) = Fcms_dot;
ydot(29) = Fcmi_dot;
ydot(30) = Fcq_dot;
ydot(31) = casub_dot;
ydot(32) = cai_dot;
ydot(33) = carel_dot;
ydot(34) = caup_dot;
ydot(35) = (1-Na_clamp)*nai_dot;
ydot(36) = 0*ki_dot;
ydot(37) = v_dot;
ydot(38) = yb_1_2_dot;
%% Adjust outputs depending on the function call

if (nargin == 3)
    output = ydot;
elseif (nargin == 4) && strcmp(runType,'ydot')
    output = ydot;
elseif (nargin == 4) && strcmp(runType,'dVm')
    currents = [v_dot];
    output = currents;
elseif (nargin == 4) && strcmp(runType,'currents')
    currents = 1/capacitance * [v_dot*capacitance ist ina_ttxs ina_ttxr...
        icat ical12 ical13 ih ik1 ikr iks ito isus...
        ibna ibca ibk inak inaca Jrel*capacitance Jup*capacitance ihna];
    output = currents;
end
% end ecc function
