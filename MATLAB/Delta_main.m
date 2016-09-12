%%Main file for loading rate inference project
%First call synthetic data function
A = [.2,.3,.2;.4,.6,.7;.4,.1,.1];
V = [0,2,4];
E_NOISE = [0,0];
E_FREQ = [5,1];
NOISE = 0;
SEQ = 50;
MEM = 5;
RES = 1;
%{
Args:
        seq_length (int): # of promoter states in sequence
        K (int): number of naive states
        res (int): number of discrete time steps per promoter state 
        w (int): memory
        A (stochastic float matrix): transition matrix
        v (float array): emission values
        v_noise (float): Scale of STD of Gaussian noise relative to 
        average state sep. Imitates inherent fluctuations in emission
        values
        exp_freq (int array): really the period of noise signals in time
        units
        exp_noise (float array): magnitude of gaussian noise signals
        relative to state separation
%}    

[fluo, fluo_interp, compound, naive] = Synthetic (SEQ, 3, RES, MEM, A, V, NOISE, E_FREQ, E_NOISE);
fluo_interp = horzcat(0,fluo_interp);
%zero negative terms
fluo_interp(fluo_interp<0) = 0;

t_interp = linspace(1,length(fluo), length(fluo_interp));

plot(1:SEQ, compound, 1:SEQ, fluo, t_interp, fluo_interp);
%%
length(fluo_interp)
%% Now take simple-minded differences approach (no filtering/smoothing)
%method yields unsatisfactory results, even when only noise is same freq
%as signal (i.e. state duration). Even worse when higher freq noise is added 
f_diffs = horzcat(0,diff(fluo_interp)).*RES;
f_2d = horzcat(0,diff(f_diffs)).*RES;
t_diffs = horzcat(0,t_interp(2:end));

%infer loading rates using diffs and mem
%l_rates = loading_rates(f_diffs, MEM, RES);
plot(1:SEQ, naive)
% histogram(l_rates)
%% Next try applying Savitzky-Golay filter before differntiating
W = 11; %width of interpolation window in t steps
HALF_W = (W-1)/2;
DEGREE = 9; %degree of poynomial to be used

%use sgolay to create filter of desired order (k) 
fluo_filt = sgolayfilt(fluo_interp,DEGREE,W);
%differentiate
smf_diffs = horzcat(0,diff(fluo_filt)).*RES;

%infer loading rates using diffs and mem
%sml_rates = loading_rates(smf_diffs, MEM, RES);

plot(1:SEQ, naive, t_interp, sml_rates)
%histogram(sml_rates)
% hold on
% histogram(naive)
