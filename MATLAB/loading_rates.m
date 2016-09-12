%Given differential of a fluorescent signal, estimate loading rates
function [l_rates] = loading_rates(fluo, w, res)
    L = length(fluo);
    %initialize l_rate vector
    W = res*w;
    l_rates = horzcat(fluo(1:W),zeros(1,L-W));
    %iterate through diff signal
    for i = (W+1):L
        u = max(l_rates(i-W),0);
        l_rates(i) = fluo(i)+u;
    end
    l_rates = l_rates;
        
    
    