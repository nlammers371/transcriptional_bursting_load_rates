
function [fluo, fluo_noise, e_agg, emissions]...
    = Synthetic (seq_length, K, res, w, A, v, v_noise, exp_freq, exp_noise)
%{
   Generates a fluorescence sequence using the given model parameters

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
        exp_freq (int array): really the period of noise signals in units state/res
        exp_noise (float array): magnitude of gaussian noise signals
        relative to state separation
    

    Returns:
        Fluorescence sequence of length seq_length, generated randomly using
        the model parameters. Gaussian noise is added at all time points.
    %}
    naive_states = transpose(zeros(seq_length,1));
    % assigns naive state 1 as the first state in the sequence
    naive_states(1) = 1;

    % generates a sequence of naive states using the transition probabilities
    for t = 2:seq_length
        naive_states(t) = randsample(1:K,1,true,A(:,naive_states(t-1)));
    end
    
    %convert to emission values
    emissions = v(naive_states);
    
    %use w to calculate aggregate emission states
    e_sum = cumsum(emissions);
    e_agg = cat(2,e_sum(1:w),(e_sum((w+1):end) - e_sum(1:(end-w))));
    
    %calculate average distance between states
    state_dist = mean(diff(v));
    %add gaussian noise term
    fluo = e_agg + normrnd(0,v_noise*state_dist,1,seq_length);
    
    %interpolate to achive desired resolution relative to promoter state
    %time scale
    t_interp = linspace(1,seq_length, res*seq_length);
    fluo_interp = interp1(1:seq_length,fluo,t_interp);
   
    %add exp noise
    e_noise = transpose(zeros(seq_length*res,1));
    for i = 1:length(exp_noise)
        err = normrnd(0,exp_noise(i)*state_dist,1,seq_length*res/exp_freq(i));
        e_noise = e_noise + repelem(err, exp_freq(i));
    end
    
    %add noise to signal
    fluo_noise = fluo_interp + e_noise;
        