function extraConvolutionDemos
    close all
    
%     step_test

    t = linspace(0,5,1000);
    Cp = AIF(t);
    
%     figure, plot(t,Cp)
    
    ti = linspace(t(1), t(end), 50);
    dti = ti(2)-ti(1);
    
    t0 = 1;
%     kEp = 5;
%     kTrans = 1/kEp;
    kEp = 5;
    kTrans = 1;
    %kTrans = 1/(1+1/kEp);
    %kTrans = 1/(1 + exp(-kEp*dti))*1/sqrt(2*dti);
    
    %e = kTrans * exp( -kEp*(t - t0) ) .* (t >= t0);
    e = kernel(t, kTrans, kEp, t0);
    figure, hold all; plot(t,e)
    
    E = kernel(t', kTrans, kEp, ti);
    Cpi = AIF(ti);
    
    figure('Position', [717   524   560   420])
    plot(repmat(t',[1 size(E,2)]), bsxfun(@times, E, Cpi)/kTrans)
    hold all
    plot(t,Cp)
    
    %figure('Position', [717   19   560   420])
    %plot(ti, sum(bsxfun(@times, E, Cpi)),'-o')
    
    figure('Position', [717   19   560   420])
    plot(ti, convolveCWithKernel(t, kTrans, kEp, ti))
    
    
%     hold all
%     plot(ti, convolveCWithKernel(ti, kTrans, kEp, ti),'-o')
%     plot(t, sum(bsxfun(@times, E, Cp)),'-')
    
    snapnow
    
    
    %% 
    %
    %
    f1 = @(x) AIF(x);
    f2 = @(x) kernel(x, kTrans, kEp, 0);
    %t = linspace(t(1),t(end),4000);
    q = quadv(@(tau) f2(tau) * f1(t-tau), t(1), t(end));
    
    figure, plot(t,q)
    hold all
    plot(ti, convolveCWithKernel(t, kTrans, kEp, ti)/300)
    
    snapnow
    
    keps = logspace(log10(0.05),log10(100),10);
    s = zeros(length(t),length(keps));
    for k = 1:length(keps)
        for i=1:length(t)
            s(i,k) = keps(k)*convolutionFromMaple(t(i), keps(k), 1, 1);
        end
    end
    figure, plot(t,s, 'LineWidth', 0.5)
    %print -depsc2 -r300 myfile.eps
    ti_old = ti;


    oversample_i = 4;
    oversample_j = 8;
    ti = linspace(ti_old(1),ti_old(end),oversample_i*length(ti_old));
    tj = linspace(ti_old(1),ti_old(end),oversample_j*length(ti_old));
    Cpj = AIF(tj);
    
    figure
    hold all
    plot(tj,Cpj*10,'LineWidth',5)
    
    for k_ep = logspace(log10(0.01),log10(1000),20);
        KTrans = k_ep;
        si = zeros(size(ti));
        for j=1:length(tj)
            si = si + Cpj(j) * KTrans * convolutionFromMapleVectorized(ti, k_ep, tj(j), oversample_j);
    %         for i=1:length(si)
    %             si(i) = si(i) + Cpi(j) * convolutionFromMaple(ti(i), 2, ti(j), oversample);
    %         end
        end
        plot(ti,si)
    end
    %figure, plot(ti,si)
    
    
    
    %convolutionDemo
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

function e = kernel(t, kTrans, kEp, t0)
    e = kTrans * exp( -kEp*bsxfun(@minus, t, t0) ) .* bsxfun(@ge, t, t0);
end

function g = convolveCWithKernel(t, kTrans, kEp, ti)
    E = kernel(t', kTrans, kEp, ti);
    Cpi = AIF(ti);
    
    g = sum(bsxfun(@times, E, Cpi));
end

