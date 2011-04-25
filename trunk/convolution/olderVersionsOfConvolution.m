

%%
function signal = convolutionForC_optimize_1(KTrans, k_ep, dt_i, Ti, dt_j, Tj, Cpi, samplingRate)
    
    % Interval length
    L = 1/samplingRate;

    % Common subexpressions
    f = k_ep*L;
    a = exp(f);
    ai= 1/a;
    b = ai - 2 + a;
    c = KTrans * samplingRate / (k_ep * k_ep);
    
    % Scale the input function (vector) for the convolution
    Ci = zeros(1,Ti);
    for i=1:Ti
        Ci(i) = c * Cpi(i);
    end

    % Compute the convolution
    signal = zeros(1,Tj);
    for i = 1:Ti
        for j = 1:Tj
            tj = dt_j * j;
            ti = dt_i * i;
            u = tj - ti;

            % Some more common terms
            g = k_ep*u;
            e = exp(-g);
            
            % Fake-branch
            if u <= -L
                s = 0;
            elseif u <= 0 
                s = e * ai - 1 + f + g; 
            elseif u <= L
                s = e * (ai - 2) + 1 + f - g;
            else
                s = e * b;
            end
            
            % Accumulate
            signal(j) = signal(j) + Ci(i) * s;
        end
    end
end


%%
function signal = convolutionForC(KTrans, k_ep, ti, tj, Cpi, samplingRate)
    
    % Interval length
    L = 1/samplingRate;
    
    Ti = numel(ti);
    Tj = numel(tj);

    % Common factors
    a = exp(k_ep*L);
    b = 1/a - 2 + a;
    c = KTrans * samplingRate / (k_ep * k_ep);
    
    % Scale the input function (vector) for the convolution
    Ci = zeros(1,Ti);
    for i=1:Ti
        Ci(i) = c * Cpi(i);
    end

    % Compute the convolution
    signal = zeros(1,Tj);
    for i = 1:Ti
        for j = 1:Tj
            u = tj(j) - ti(i);
            if u <= -L
                s = 0;
            elseif u <= 0 
                s = exp(-k_ep*(L + u)) - 1 + k_ep*(u + L); 
            elseif u <= L
                s = exp(-k_ep*(L + u)) - 2*exp(-k_ep*u) + 1 + k_ep*(L - u);
            else
                s = exp(-k_ep*u) * b;
            end
            signal(j) = signal(j) + Ci(i) * s;
        end
    end
end


%%
function signal = convolutionInnerLoop(t, k_ep, t_0, samplingRate)
    L = 1/samplingRate;
    T = length(t);
    
    a = exp(k_ep*L);
    b = 1/a - 2 + a;
    
    signal = zeros(1,T);
    for j = 1:T
        u = t(j) - t_0;
        if u <= -L
            s = 0;
        elseif u <= 0 
            s = exp(-k_ep*(L + u)) - 1 + k_ep*(u + L); 
        elseif u <= L
            s = exp(-k_ep*(L + u)) - 2*exp(-k_ep*u) + 1 + k_ep*(L - u);
        else
            s = exp(-k_ep*u) * b;
        end
        s = s * samplingRate / (k_ep * k_ep);
        
        signal(j) = s;
    end
end


%%
function s = convolutionFromMapleVectorized(t, k, t_0, oversamplingFactor)
    x = t - t_0;
    L = 1/oversamplingFactor;
    s = zeros(size(t));
    
    ind_1 = (x > -L  &  x <= 0);
    ind_2 = (x >  0  &  x <= L);
    ind_3 = (x > L);
    
    a = exp(k*L);
    s(ind_1) = exp(-k*(L + x(ind_1))) - 1 + k*(x(ind_1) + L); 
    s(ind_2) = exp(-k*(L + x(ind_2))) - 2*exp(-k*x(ind_2)) + 1 + k*(L - x(ind_2));
    s(ind_3) = exp(-k*x(ind_3)) * (1/a - 2 + a);

    s = s * oversamplingFactor / (k * k);
end

function s = convolutionFromMaple(t, k, t_0, oversamplingFactor)
    x = t - t_0;
    L = 1/oversamplingFactor;
    if x <= -L
        s = 0;
    else
        if x > L
            %s = exp(-k*(L + x)) - 2*exp(-k*x) + exp(-k*(-L + x));
            %s = exp(-k*x) * (exp(-k*L) - 2 + exp(k*L));
            a = exp(k*L);
            s = exp(-k*x) * (1/a - 2 + a);
            %s = exp(-k*x) * sinh(k*L/2)^2;
        elseif x <= 0
            s = exp(-k*(L + x)) - 1 + k*(x + L); 
        elseif x <= L
            s = exp(-k*(L + x)) - 2*exp(-k*x) + 1 + k*(L - x);
        end
        s = s * oversamplingFactor / (k * k);
    end
end

function s = convolutionFromMapleOld(t, k, n)
    if t <= -1 + n
        s = 0;
    elseif t <= n
        s = (exp(((-1 - t + n) * k)) - 1 + (k * (t - n + 1))) / (k ^ 2); 
    elseif t <= 1 + n
        s = (exp(((-1 - t + n) * k)) - 2 * exp(-(k * (t - n))) + 1 + (k * (1 - t + n))) / (k ^ 2);
    elseif 1 + n < t
        s = (exp(-(k * (t - n - 1))) + exp(((-1 - t + n) * k)) - 2 * exp(-(k * (t - n)))) / (k ^ 2);
    end

%     cg7return = piecewise(
%         t <= -1 + n, 0, 
%         t <= n, (exp(((-1 - t + n) * k)) - 0.1e1 + (k * (t - n + 1))) / (k ^ 2), 
%         t <= 1 + n, (exp(((-1 - t + n) * k)) - 0.2e1 * exp(-(k * (t - n))) + 0.1e1 + (k * (1 - t + n))) / (k ^ 2), 
%         1 + n < t, (exp(-(k * (t - n - 1))) + exp(((-1 - t + n) * k)) - 0.2e1 * exp(-(k * (t - n)))) / (k ^ 2)
%     );
end

function e = kernelScalarVersion(t, kTrans, kEp, t0)
    e = kTrans * exp( -kEp*(t - t0) ) .* (t >= t0);
end

function e = kernelAltVersion(t, kTrans, kEp)
    e = kTrans * exp( -kEp*(t) );
end

function e = kernel(t, kTrans, kEp, t0)
    e = kTrans * exp( -kEp*bsxfun(@minus, t, t0) ) .* bsxfun(@ge, t, t0);
end

function g = convolveCWithKernel(t, kTrans, kEp, ti)
    E = kernel(t', kTrans, kEp, ti);
    Cpi = AIF(ti);
    
    g = sum(bsxfun(@times, E, Cpi));
end
