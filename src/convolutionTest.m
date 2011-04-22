%% Convolution Function
% Convolution of C with decaying exponential kernel.
%
% $$\int_0^t K^{trans}e^{-k_{ep} \tau} C(t-\tau - t_0) d\tau$$
% 

%% Syntax
%
% * *prototypeConvolution* -- Calling with no arguments runs a demo.
% * g = *prototypeConvolution*(kTrans, kEp, t0, t) -- Typical call.
%       
% 
%% Description
% Implements the convolution of C with decaying exponential kernel, as described in
%       Suyash et al.  The output of this function is a time-series representation of the 
%       pixel intensity over time.  It is also the solution to the first order ODE
%       modeling the flow of tracer from the bloodstream to the tissue and back.
%       
%       Int( C(tau) * exp(t - tau) dtau, 0, t)
% 
% $$e^{\pi i} + 1 = 0$$ 
% 
%   Foo  bar
% 
% $$e^{\pi i} + 1 = 0$$
% 
% 
% # ITEM1
% # ITEM2
% 
% $$\int_0^t K^{trans}e^{-k_{ep} \tau} C(t-\tau) d\tau$$
% 

%%
% 
%
%   Arguments:
%       kTrans      Scale parameter (beta_1)
%       kEp         Shape parameter (1/beta_2)
%       t0          Shift or Delay parameter (beta_3)
%       t           Temporal coordinate vector or Time-points (t)
%
%   Output:
%       g           Pixel intensity time series (g(x,t) in TMI paper)
%

%% Arguments
%  kTrans      Scale parameter (beta_1)
% kEp         Shape parameter (1/beta_2)
% t0          Shift or Delay parameter (beta_3)
% t           Temporal coordinate vector or Time-points (t)
%
%%  Output
% g           Pixel intensity time series (g(x,t) in TMI paper)
%

function prototypeConvolution(~)
    
    if ~nargin, demo, return, end

end


%% Demo
%
function demo
    
    close all

    t = linspace(0,5,1000);
    Cp = breastCp(t);
    
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
    hold all; plot(t,e)
    
    E = kernel(t', kTrans, kEp, ti);
    Cpi = breastCp(ti);
    
%     figure('Position', [717   524   560   420])
%     plot(repmat(t',[1 size(E,2)]), bsxfun(@times, E, Cpi)/kTrans)
%     hold all
%     plot(t,Cp)
%     
%     %figure('Position', [717   19   560   420])
%     %plot(ti, sum(bsxfun(@times, E, Cpi)),'-o')
    
%     figure('Position', [717   19   560   420])
%     plot(ti, convolveCWithKernel(t, kTrans, kEp, ti))
    
    
    %hold all
    %plot(ti, convolveCWithKernel(ti, kTrans, kEp, ti),'-o')
    %plot(t, sum(bsxfun(@times, E, Cp)),'-')
    
%     snapnow
    
    
    %% 
    %
    %
    f1 = @(x) breastCp(x);
    f2 = @(x) kernel(x, kTrans, kEp, 0);
    %t = linspace(t(1),t(end),4000);
    q = quadv(@(tau) f2(tau) * f1(t-tau), t(1), t(end));
    
    figure, plot(t,q)
    hold all
    plot(ti, convolveCWithKernel(t, kTrans, kEp, ti)/300)
    
    snapnow
    
    s = zeros(size(t));
    for i=1:length(s)
        s(i) = convolutionFromMaple(t(i), 5, 2, 1);
    end
    figure, plot(t,s)
    ti_old = ti;


    oversample_i = 2;
    oversample_j = 8;
    ti = linspace(ti_old(1),ti_old(end),oversample_i*length(ti_old));
    tj = linspace(ti_old(1),ti_old(end),oversample_j*length(ti_old));
    Cpj = breastCp(tj);
    
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
    
    
    
    convolutionOuterLoop

    '';
end


%%
function signal = convolutionOuterLoop()
    t0 = 0; 
    tf = 5;
    T = 50;

    oversample_i = 8;
    oversample_j = 2;
    
    Ti = oversample_i*T;
    Tj = oversample_j*T;
    
    ti = linspace(t0, tf, Ti);
    tj = linspace(t0, tf, Tj);

    Cpi = breastCp(ti);
    
    dt_i = ti(2) - ti(1);
    dt_j = tj(2) - tj(1);

    figure
    hold all
    plot(ti, Cpi*10, 'LineWidth', 5)
    
    for k_ep = logspace(log10(0.01), log10(1000), 20);
        KTrans = k_ep;
%         
%         sj = zeros(1,Tj);
%         for i=1:Ti
%             signal_part = convolutionInnerLoop(tj, k_ep, ti(i), oversample_i);
%             sj = sj + Cpi(i) * KTrans * signal_part;
%         end

        signal = dce_mri_mex(KTrans, k_ep, dt_i, Ti, dt_j, Tj, Cpi, oversample_i);
        plot(tj, signal)
    end
end


%%
%
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
    Cpi = breastCp(ti);
    
    g = sum(bsxfun(@times, E, Cpi));
end

%% Population model for Breast Input function
function Cp = breastCp(t,dt)
    if (nargin < 2) dt = 32/60; end;

    t = 60*t;
    dt = 60*dt;

    alpha = 4.0;

    %Cp = Cps([0.10 dt+48 alpha 20.0 1000.0],t); 

    Cp = Cpg([1.0 dt alpha 5.0],t) + ...
         Cpg([0.2 dt+28 alpha 10.0],t) + ...
         Cpg([0.10 dt+38 alpha 20.0],t) + ...
         Cps([0.10 dt+48 alpha 20.0 1000.0],t);

    Cp = 5.0*Cp;

end

% param(1) = Cp0
% param(2) = Delta
% param(3) = alpha
% param(4) = tau
%
function f = Cpg(param,t)

    Cp0 = abs(param(1));
    Delta = param(2);
    alpha = max(abs(param(3)),1);
    tau = abs(param(4));

    tprime = t-Delta;

    alphaprime = alpha-1;
    norm = exp(-alphaprime)*(tau*alphaprime)^alphaprime;

    if (length(t) == 1)
        if (t < Delta) 
            f = 0;
        else
            f = (Cp0/norm)*(tprime^alphaprime)*exp(-tprime/tau);
        end;
    else
        idx = find(t>=Delta);
        f = zeros(size(t));
        f(idx) = (Cp0/norm)*(tprime(idx).^alphaprime).*exp(-tprime(idx)/tau);
    end;

end

% param(1) = Cp0
% param(2) = delta
% param(3) = alpha
% param(4) = tau
% param(5) = T
%
function f = Cps(param,t)

    Cp0 = abs(param(1));
    Delta = param(2);
    alpha = max(abs(param(3)),1);
    tau = abs(param(4));
    T = abs(param(5));

    tprime = t-Delta;

    norm = gamma(alpha);

    xi = 1/tau-1/T;
    gam = norm;

    if (length(t) == 1)
        if (t < Delta) 
            f = 0;
        else
            f = (Cp0/norm)*gam*exp(-tprime/T)*gammaincc(xi*tprime,alpha);
        end;
    else
        idx = find(t>=Delta);
        f = zeros(size(t));
        f(idx) = (Cp0/norm)*gam*exp(-tprime(idx)/T).*gammaincc(xi*tprime(idx),alpha);
    end;

end

%% Corrected gammainc?
function y = gammaincc(x, a)
    y = gammainc(x, a);
end
