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

function step_test
    X = 100;
    Y = 50;
    
    C = zeros(X,Y);
    
    x0 =   0; y0 =  0;
    xf = X-1; yf = 35;
    dx = xf - x0;
    dy = yf - y0;
    
    %xt = 0; yt = 0;
    xt = floor(dx/2); yt = 0;%floor(dy/2);
    %xt = 0; yt = floor(dy/2);
    y = y0;
    for x = x0:xf
        C(x+1,y+1) = 1;
        
        yt = yt + dy;
        if yt > xt
            xt = xt + dx;
            y = y + 1;
        end
    end
    
    figure
    imagesc(C','XData',[0 X-1]-0/2,'YData',[0 Y-1]-0/2), axis xy
    
    set(gca,'XTick', (0:X-1)-1/2)
    set(gca,'YTick', (0:Y-1)-1/2)
    set(gca,'XGrid','on')
    set(gca,'YGrid','on')
    set(gca,'XColor',[0.5,0.6,0.8])
    set(gca,'YColor',[0.5,0.6,0.8])
    set(gca,'GridLineStyle', '-')    
    '';
end

%% Demo
%
function demo
%     
%     close all
%     
%     step_test
% 
%     t = linspace(0,5,1000);
%     Cp = AIF(t);
%     
% %     figure, plot(t,Cp)
%     
%     ti = linspace(t(1), t(end), 50);
%     dti = ti(2)-ti(1);
%     
%     t0 = 1;
% %     kEp = 5;
% %     kTrans = 1/kEp;
%     kEp = 5;
%     kTrans = 1;
%     %kTrans = 1/(1+1/kEp);
%     %kTrans = 1/(1 + exp(-kEp*dti))*1/sqrt(2*dti);
%     
%     %e = kTrans * exp( -kEp*(t - t0) ) .* (t >= t0);
%     e = kernel(t, kTrans, kEp, t0);
%     figure, hold all; plot(t,e)
%     
%     E = kernel(t', kTrans, kEp, ti);
%     Cpi = AIF(ti);
%     
%     figure('Position', [717   524   560   420])
%     plot(repmat(t',[1 size(E,2)]), bsxfun(@times, E, Cpi)/kTrans)
%     hold all
%     plot(t,Cp)
%     
%     %figure('Position', [717   19   560   420])
%     %plot(ti, sum(bsxfun(@times, E, Cpi)),'-o')
%     
%     figure('Position', [717   19   560   420])
%     plot(ti, convolveCWithKernel(t, kTrans, kEp, ti))
%     
%     
% %     hold all
% %     plot(ti, convolveCWithKernel(ti, kTrans, kEp, ti),'-o')
% %     plot(t, sum(bsxfun(@times, E, Cp)),'-')
%     
%     snapnow
%     
%     
%     %% 
%     %
%     %
%     f1 = @(x) AIF(x);
%     f2 = @(x) kernel(x, kTrans, kEp, 0);
%     %t = linspace(t(1),t(end),4000);
%     q = quadv(@(tau) f2(tau) * f1(t-tau), t(1), t(end));
%     
%     figure, plot(t,q)
%     hold all
%     plot(ti, convolveCWithKernel(t, kTrans, kEp, ti)/300)
%     
%     snapnow
%     
%     keps = logspace(log10(0.05),log10(100),10);
%     s = zeros(length(t),length(keps));
%     for k = 1:length(keps)
%         for i=1:length(t)
%             s(i,k) = keps(k)*convolutionFromMaple(t(i), keps(k), 1, 1);
%         end
%     end
%     figure, plot(t,s, 'LineWidth', 0.5)
%     %print -depsc2 -r300 myfile.eps
%     ti_old = ti;
% 
% 
%     oversample_i = 4;
%     oversample_j = 8;
%     ti = linspace(ti_old(1),ti_old(end),oversample_i*length(ti_old));
%     tj = linspace(ti_old(1),ti_old(end),oversample_j*length(ti_old));
%     Cpj = AIF(tj);
%     
%     figure
%     hold all
%     plot(tj,Cpj*10,'LineWidth',5)
%     
%     for k_ep = logspace(log10(0.01),log10(1000),20);
%         KTrans = k_ep;
%         si = zeros(size(ti));
%         for j=1:length(tj)
%             si = si + Cpj(j) * KTrans * convolutionFromMapleVectorized(ti, k_ep, tj(j), oversample_j);
%     %         for i=1:length(si)
%     %             si(i) = si(i) + Cpi(j) * convolutionFromMaple(ti(i), 2, ti(j), oversample);
%     %         end
%         end
%         plot(ti,si)
%     end
%     %figure, plot(ti,si)
%     
%     
    step_test
    convolutionDemo

    '';
end


%%
function convolutionDemo()
    t0 = 0; 
    tf = 5;
    T = 50;

    oversample_i = 4;
    oversample_j = 8;
    
    Ti = oversample_i*T;
    Tj = oversample_j*T;
    
    dt_i = (tf - t0) / (Ti);
    dt_j = (tf - t0) / (Tj);
    
    ti = (0:Ti-1) * dt_i;
    tj = (0:Tj-1) * dt_j;
    
%     ti = linspace(t0, tf, Ti);
%     tj = linspace(t0, tf, Tj);

    Cpi = AIF(ti);
    
%     dt_i = ti(2) - ti(1);
%     dt_j = tj(2) - tj(1);

    figure
    hold all
    plot(ti, Cpi*10, 'LineWidth', 5)
    
    %for k_ep = logspace(log10(0.01), log10(1000), 20);
    for k_ep = logspace(log10(0.01), log10(1000), 3);
        KTrans = k_ep;
%         
%         sj = zeros(1,Tj);
%         for i=1:Ti
%             signal_part = convolutionInnerLoop(tj, k_ep, ti(i), oversample_i);
%             sj = sj + Cpi(i) * KTrans * signal_part;
%         end

        %signal = convolutionForC(KTrans, k_ep, ti, tj, Cpi, oversample_i);
        %signal = convolutionForC_optimize_2(KTrans, k_ep, dt_i, Ti, dt_j, Tj, Cpi, oversample_i);
        signal = convolutionForC_optimize_2(KTrans, k_ep, dt_i, Ti, dt_j, Tj, Cpi, oversample_i);
        plot(tj, signal)
    end
end
%%

%%%%%%%%%%%%%%%%%%  C    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signal = convolutionForC_optimize_3(KTrans, k_ep, dt_i, Ti, dt_j, Tj, Cpi, samplingRate)
    
    % Interval length
    L = 1/samplingRate;

    % Common subexpressions
    f = k_ep*L;
    a = exp(f);
    ai= 1/a;
    b = ai - 2 + a;
    c = KTrans * samplingRate / (k_ep * k_ep);
    
    % Scale the input function (vector) for the convolution
%     Ci = zeros(1,Ti);
%     for i=1:Ti
%         Ci(i) = c * Cpi(i);
%     end
    dh_i = dt_i/dt_j;
    %fn   = floor(-L/dt_j);
    %fp   = floor( L/dt_j);     % Not same as -floor(L/dt_j)

    % Do not use round() in cuda. It maps to an 8-instruction sequence.
    % Use rint() instead since it maps to a 1-instruction sequence.
    
    L_dt    = L / dt_j;
    L_dt_r  = round(L_dt * 1e9) * 1e-9;  % In case L is a multiple of dt_j

    j00 = 0;
    j01 = ceil(-L_dt_r);
    j02 = 0;
    j03 = ceil( L_dt_r);
    j04 = Tj;

%     x0 = -dx;
%     y0 = -dy;
    
    [dy, dx] = rat(L_dt, 1e-9);
    x0 = 0; y0 = 0;
    x = x0; y = y0;
    

    % dy < dx if L_dt < 1
    x = x + dx;
    y = y + dy;

    error 'implement me'
    
    
    z0 = -L_dt_r;
    dz = dt_i;
    
    z = z0;
    z = z + dz;
    
    
    
    
    % Compute the convolution
    signal = zeros(1,Tj);
    hi = 0 + Ti*eps(dh_i);
    ti = 0;
    for i = 0:Ti-1
        Ci = c * Cpi(i+1);
        
        ti2 = dt_i * i;
        %ti


%         s1 = 0;
%         s2 = e * ai - 1 + f + g;
%         s3 = e * (ai - 2) + 1 + f - g;
%         s4 = e * b;

%         j0 = 0;
%         j1 = floor(hi - L_dt);
%         j2 = floor(hi);
%         j3 = floor(hi + L_dt);
%         j4 = Tj;
        
        j0 = 0;
        j1 = ceil((ti - L)/dt_j);
        j2 = ceil((ti    )/dt_j);
        j3 = ceil((ti + L)/dt_j);
        j4 = Tj;
        
        
        
%         j1 = max(j0, min(j4, j1)); 
%         j2 = max(j0, min(j4, j2)); 
%         j3 = max(j0, min(j4, j3)); 
        
        % Branch 2
        % s2 = e * ai - 1 + f + g;
        tj  = dt_j * j1;
        u   = tj - ti;
        g   = k_ep*u;
        dg  = k_ep*dt_j;
        dgC = dg*Ci;
        e   = exp(-g);
        me  = exp(-dg);
        
        c0 = (-1 + f)*Ci + g*Ci;
        c1 = ai * Ci;
        for j = j1:j2-1
            %e = exp(-g);
            %if j > 0, signal(j+1) = signal(j+1) + c1*e + c0; end
            %g = g + dg;
            
            if j >= 0 && u > -L
                tj = dt_j * j;
                u = tj - ti;
                g = k_ep*u;
                E = exp(-g);
                s2 = (E * ai - 1 + f + g)*Ci;
                
                %c0 = (-1 + f)*Ci + g*Ci;
                
                %s2 = Ci * (exp(-g-f) - 1 + f + g);
                s2_alt = c1*e + c0;
                
                signal(j+1) = signal(j+1) + s2;
                
                %if abs(s2 - s2_alt) > 1e2
                %if abs(s2 - s2_alt) > 1e0
                if abs(s2 - s2_alt) > 1e-2  ||  u < -L
                    % tj - ti == u >= -L
                    % tj == ti + u >= ti - L
                    % j*dt_j == ti + u >= ti - L
                    % j == (ti + u)/dt_j >= (ti - L)/dt_j
                    % --> j1 == ceil(ti - L)/dt_j
                    '';
                end
            end
            c0 = c0 + dgC;
            e = e * me;
        end
        
        % Branch 3
        % s3 = e * (ai - 2) + 1 + f - g;
        tj  = dt_j * j2;
        u   = tj - ti;
        g   = k_ep*u;
        %dg  = k_ep*dt_j;
        e   = exp(-g);
        %me  = exp(-dg);
        %c0 = (-1 + f)*Ci + g*Ci;
        %c1 = ai * Ci;
        
        % s3 = e * (ai - 2) + 1 + f - g;
        c0  = (1 + f - g) * Ci;
        c1  = (ai - 2) * Ci;
        for j = j2:j3-1
            if j < length(signal)
                signal(j+1) = signal(j+1) + c1*e + c0;
            end
            c0 = c0 - dgC;
            e = e * me;
        end

        
        % Branch 4
        % s4 = e * b;
        %tj  = dt_j * j2;
%         tj  = dt_j * j3;
%         u   = tj - ti;
%         g   = k_ep*u;
        %dg  = k_ep*dt_j;
        %e   = exp(-g);
%         e   = exp(-g);
%         me  = exp(-dg);
        
%         sj3 = signal(min(end,j3+1));


        c0  = 0;
        c1  = b * Ci;
        for j = j3:j4-1
            if j < length(signal)
                signal(j+1) = signal(j+1) + c1*e + c0;
            end
            e = e * me;
        end
        
%         j3_float = (ti + L)/dt_j;
%         dj = j3 - j3_float;
%         signal(min(end,j3+1)) = sj3 + c1*((0)*exp(-g) + (1-dj)*exp(-g-dg));
        
        
        ti - ti2
        ti = ti + dt_i;
%         hi = hi + dh_i;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function signal = convolutionForC_optimize_2(KTrans, k_ep, dt_i, Ti, dt_j, Tj, Cpi, samplingRate)
    
    % Interval length
    L = 1/samplingRate;

    % Common subexpressions
    f = k_ep*L;
    a = exp(f);
    ai= 1/a;
    b = ai - 2 + a;
    c = KTrans * samplingRate / (k_ep * k_ep);
    
    % Scale the input function (vector) for the convolution
%     Ci = zeros(1,Ti);
%     for i=1:Ti
%         Ci(i) = c * Cpi(i);
%     end
    dh_i = dt_i/dt_j;
    %fn   = floor(-L/dt_j);
    %fp   = floor( L/dt_j);     % Not same as -floor(L/dt_j)

    L_dt   = L/dt_j;
    
    
    % Compute the convolution
    signal = zeros(1,Tj);
    hi = 0 + Ti*eps(dh_i);
    ti = 0;
    for i = 0:Ti-1
        Ci = c * Cpi(i+1);
        %ti = dt_i * i;
        
        % Shouldn't these increment at the end of the for loop?
        ti = ti + dt_i;
        hi = hi + dh_i;
        
        %-L/dt_j
        
%         j1 = floor(ti/dt_j) + floor(-L/dt_j);
%         j1 = floor(i*dt_i/dt_j) + floor(-L/dt_j);
        
        
%         j0 = 0;
%         j1 = floor((ti - L)/dt_j);
%         j2 = floor((ti    )/dt_j);
%         j3 = floor((ti + L)/dt_j);
%         j4 = Tj;

        j0 = 0;
        j1 = floor(hi - L_dt);
        j2 = floor(hi);
        j3 = floor(hi + L_dt);
        j4 = Tj;
        
        
        
        for j = 0:Tj-1
            tj = dt_j * j;
            u = tj - ti;

            % Some more common terms
            g = k_ep*u;
            e = exp(-g);
            
            % u <= -L
            % tj - ti <= -L
            %  j*dt_j <= -L + ti
            %       j <= (-L + ti)/dt_j;

            % u <= -L   -->     j <= (ti - L)/dt_j;
            % u <= 0    -->     j <= (ti    )/dt_j;
            % u <= L    -->     j <= (ti + L)/dt_j;
            branch = 0;
            if j <= j1
                branch = 1;
            elseif j <= j2
                branch = 2;
            elseif j <= j3
                branch = 3;
            elseif j <= j4
                branch = 4;
            else
                branch = 5;
            end
                
            s1 = 0;
            s2 = e * ai - 1 + f + g;
            s3 = e * (ai - 2) + 1 + f - g;
            s4 = e * b;
            
            % Fake-branch
            if u < -L
                s = s1;
                otherBranch = 1;
            elseif u < 0 
                s = s2; 
                otherBranch = 2;
            elseif u < L
                s = s3;
                otherBranch = 3;
            else
                s = s4;
                otherBranch = 4;
            end

            if branch ~= otherBranch
                fprintf(['Warning: expected branch %d, was really in %d: ' ...
                'delta=%d, j = %d, (%d,%d,%d,%d,%d), u=%g, L=%g\n'], ...
                otherBranch, branch,branch-otherBranch,j,j0,j1,j2,j3,j4,u,L);
                si = [s1, s2, s3, s4];%./s;
                si(otherBranch)-si(branch)
            '';
            end
            
            % Accumulate
            signal(j+1) = signal(j+1) + Ci * s;
        end
    end
end


% {
% function s = convolutionFromMapleVectorized(t, k, t_0, oversamplingFactor)
%     x = t - t_0;
%     L = 1/oversamplingFactor;
%     s = zeros(size(t));
%     
%     ind_1 = (x > -L  &  x <= 0);
%     ind_2 = (x >  0  &  x <= L);
%     ind_3 = (x > L);
%     
%     a = exp(k*L);
%     s(ind_1) = exp(-k*(L + x(ind_1))) - 1 + k*(x(ind_1) + L); 
%     s(ind_2) = exp(-k*(L + x(ind_2))) - 2*exp(-k*x(ind_2)) + 1 + k*(L - x(ind_2));
%     s(ind_3) = exp(-k*x(ind_3)) * (1/a - 2 + a);
% 
%     s = s * oversamplingFactor / (k * k);
% end
% 
% function s = convolutionFromMaple(t, k, t_0, oversamplingFactor)
%     x = t - t_0;
%     L = 1/oversamplingFactor;
%     if x <= -L
%         s = 0;
%     else
%         if x > L
%             %s = exp(-k*(L + x)) - 2*exp(-k*x) + exp(-k*(-L + x));
%             %s = exp(-k*x) * (exp(-k*L) - 2 + exp(k*L));
%             a = exp(k*L);
%             s = exp(-k*x) * (1/a - 2 + a);
%             %s = exp(-k*x) * sinh(k*L/2)^2;
%         elseif x <= 0
%             s = exp(-k*(L + x)) - 1 + k*(x + L); 
%         elseif x <= L
%             s = exp(-k*(L + x)) - 2*exp(-k*x) + 1 + k*(L - x);
%         end
%         s = s * oversamplingFactor / (k * k);
%     end
% end
% 
% function e = kernel(t, kTrans, kEp, t0)
%     e = kTrans * exp( -kEp*bsxfun(@minus, t, t0) ) .* bsxfun(@ge, t, t0);
% end
% 
% function g = convolveCWithKernel(t, kTrans, kEp, ti)
%     E = kernel(t', kTrans, kEp, ti);
%     Cpi = AIF(ti);
%     
%     g = sum(bsxfun(@times, E, Cpi));
% end
% }

