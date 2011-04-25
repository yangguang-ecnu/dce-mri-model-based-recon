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

%% Shortened Demo
%
function shortenedDemo
    %convolutionOuterLoop;
    convolutionOuterLoop8x4;
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
    
    
    % Moved to shortenedDemo
    %convolutionOuterLoop;
    %convolutionOuterLoop4x8;
    
    shortenedDemo;

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

        signal = dce_mri_mex(single(KTrans), single(k_ep), single(dt_i), int32(Ti), single(dt_j), int32(Tj), single(Cpi), single(oversample_i));
        % plot(tj, signal)
        plot(tj, real(signal(:)))
    end
end


%%
function signal = convolutionOuterLoop8x4()
    % Image dimensions
%     X = 26;
%     Y = 26;
    X = 8;
    Y = 8;

    % Initialize KTrans and k_ep to some arbitrary values
    k_ep    = zeros(X,Y,'single');
    KTrans  = zeros(X,Y,'single');
    for x = 0:X-1
        for y = 0:Y-1
            u = x/(X + eps);
            v = y/(Y + eps);
            
            w = 1/((1 - u)*(1 - v))^2;
            w = w * 1e-1;
            
            k_ep(x+1,y+1)   = w;
            KTrans(x+1,y+1) = w/10;
        end
    end
    
    [mX,mY] = meshgrid(linspace(-3,3,X),linspace(-3,3,Y)); Z = peaks(mX,mY); 
    %Z = Z.^2;
    %Z = Z - min(Z(:)) + 1e-6;
    %Z = Z - getQuantile(Z(:), 0.05);
    %Z = max(Z, 1e-4);
    
    b = 1; c = 4; f = @(x) 1/b * x./(1 + (x*c-c).^2); %plot(x, f(x))
    Z = f(Z);
    Z = Z * 10;
    Z = max(Z, 1e-1);
    
    %figure, imagesc(Z'), axis image; colorbar
    
    k_ep = Z;
    KTrans = Z/10;

    
    figure, imagesc(Z), axis image; colorbar, title('k_ep')
    
    
    % Init time variables
    t0 = 0; 
    tf = 5;
    T  = 50;

    oversample_i = 8;
    oversample_j = 2;
    
    Ti = oversample_i*T;
    Tj = oversample_j*T;
    
    dt_i = (tf - t0) / Ti;
    dt_j = (tf - t0) / Tj;
    
    %ti = linspace(t0, tf, Ti);
    %tj = linspace(t0, tf, Tj);
    ti = (0:Ti-1) * dt_i;
    tj = (0:Tj-1) * dt_j;
%     tj = linspace(t0, tf, Tj);
    
    %dt_i = ti(2) - ti(1);
    %dt_j = tj(2) - tj(1);

    % Init Cpi vector
    Cpi = breastCp(ti);
    
    
    % Convert to 32-bit floats for everything inputted into dce_mri_mex()
    KTrans  = single(KTrans);
    k_ep    = single(k_ep);
    dt_i    = single(dt_i);
    Ti      = int32(Ti);
    dt_j    = single(dt_j);
    Tj      = int32(Tj); 
    Cpi     = single(Cpi);
    oversample_i = single(oversample_i);
    
%     ti = single(ti);
%     tj = single(tj);
    
    fprintf('Beginning the timings\n')

    % Reference implementation
    signal_0 = zeros(X,Y,Tj,'single');
	tj_vec = single(0:Tj-1) * dt_j; 
    tj_vec = permute(tj_vec, [1 3 2]);  % Make it 1x1xTj
    tic
    for x = 1:X
        for y = 1:Y
            sj = zeros(1,1,Tj,'single');
            for i = 1:Ti
                %si = si + Cpj(j) * KTrans * convolutionFromMapleVectorized(ti, k_ep, tj(j), oversample_j);
                sj = sj + Cpi(i) * KTrans(x,y) * convolutionFromMapleVectorized(tj_vec, k_ep(x,y), ti(i), oversample_i);
                %signal_0(x,y,:) = dce_mri_mex( KTrans(x,y), k_ep(x,y), dt_i, Ti, dt_j, Tj, Cpi, oversample_i );
            end
            signal_0(x,y,:) = sj;
        end
    end
    fprintf('Reference implementation (host): '), toc

    
    % Scalar version of mex-wrapper call
    signal_1 = zeros(X,Y,Tj,'single');
    tic
    for x = 1:X
        for y = 1:Y
            signal_1(x,y,:) = dce_mri_mex( KTrans(x,y), k_ep(x,y), dt_i, Ti, dt_j, Tj, Cpi, oversample_i );
        end
    end
    fprintf('Scalar version of mex-wrapper call: '), toc

    
    
    
    
    % Matrix version of mex-wrapper call
    tic
    signal_2 = dce_mri_mex(KTrans, k_ep, dt_i, Ti, dt_j, Tj, Cpi, oversample_i);
    fprintf('Matrix version of mex-wrapper call: '), toc
    %signal_2 = signal_1 + randn(X,Y,Tj)*1e-3;  % Fake output

    
    signal_1 = real(signal_1);
    signal_2 = real(signal_2);
    
    signal = signal_2;  % Final return value of this function

    
    % Compare matrix versus scalar calls
    RMSE  = norm(signal_1(:) - signal_2(:)) / sqrt(X*Y*double(Tj));
    nRMSE = RMSE / norm(signal_1(:));
    
    fprintf('Matrix versus scalar:  RMSE: %g\n',  RMSE)
    fprintf('Matrix versus scalar: nRMSE: %g\n', nRMSE)
    

    % Compare cuda versus matlab
    RMSE  = norm(signal_2(:) - signal_0(:)) / sqrt(X*Y*double(Tj));
    nRMSE = RMSE / norm(signal_0(:));
    
    fprintf('Cuda versus matlab:  RMSE: %g\n',  RMSE)
    fprintf('Cuda versus matlab: nRMSE: %g\n', nRMSE)
    

    % Show curves

    figure, hold all; plot(ti, Cpi, 'LineWidth', 5)
    signal_transposed = permute(signal_0, [3 1 2]);
    plot(tj, signal_transposed(:,:))
    title('Reference curves using matlab prototype (CPU)')
        
    figure, hold all; plot(ti, Cpi, 'LineWidth', 5)   
    signal_transposed = permute(signal_2, [3 1 2]);
    plot(tj, signal_transposed(:,:))
    title('Matrix version kernel code output (GPU)')

    figure, hold all; plot(ti, 0*Cpi)
    signal_transposed = permute(signal_2 - signal_0, [3 1 2]);
    plot(tj, signal_transposed(:,:))
    title('Error difference time curves (GPU vs CPU)')
    
    
    % Show reconstructed time series of images (movie)
    fig = figure;
%     lgsignal = log10(abs(signal));
%     crange = [-2, max(lgsignal(:))];
    lgsignal = ((signal));
    crange = [0, max(lgsignal(:))];
    for t = 1:Tj
        if ishandle(fig)
            set(0, 'CurrentFigure', fig)
        else
            break
        end
        imagesc(lgsignal(:,:,t)', crange)
        axis image xy
        colorbar
        title(['t = ' num2str(t)])
        drawnow
        pause(1/60)
    end 
    
%     [mX,mY] = meshgrid(1:X,1:Y);
%     Z = peaks(mX,mY);

    
	aspectRatio = 1.33/1;
    T0 = 10;     Tf = ceil(single(Tj)/2);

    nT = Tf - T0 + 1;
    % Tx / Ty <= 2.5
    %   Ty >= Tx / 2.5
    % Tx * Ty >= nT
    %   Tx >= nT / Ty
    % Ty >= Tx / 2.5 >= nT / (2.5 * Ty)
    %   Ty^2 >= nT / 2.5
    %   Ty >= sqrt(nT/2.5)
    %   Tx >= nT / Ty
    
    
    %nT = 293; aspectRatio = 3/1; Ty = round(sqrt(nT/aspectRatio)); Tx = round((nT/Ty)); [Tx Ty Tx*Ty Ty*aspectRatio Tx/Ty]
    Ty = round(sqrt(nT/aspectRatio)); 
    Tx = round((nT/Ty)); 
    Tf = T0 + Tx*Ty - 1;
    %[Tx Ty Tx*Ty Ty*aspectRatio Tx/Ty]

    sampling = double(signal(:,:,T0:Tf));
    a = 0.05;
    p = [a/2, 1-a/2];
    range = getQuantile(sampling(:), p);


    %figure, plot(linspace(0,1,ns),s, p, range, 'o', 'MarkerSize', 14)
    
    sampling = single(mat2gray(sampling, range));
    figure, montage(permute(sampling,[2 1 4 3]), 'size', [Ty Tx]), colormap(jet)
    figure, montage(permute(sampling,[2 1 4 3]), 'size', [Ty Tx]), colormap(bone)
    
    %figure, imagesc(k_ep'), axis image; colorbar
end

%%
function q = getQuantile(x, p)
    s = sort(x(:,:));
    ns = size(s,1);
    ind = min(ns, 1+max(0, floor(ns * p(:)')));
    q = s(ind,:);
    
    % getQuantile([sampling(:), 2*sampling(:)],[a/2 0.25 0.5 0.75 1-a/2])
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
    %s = exp(-k*x) * (1/a - 2 + a)
    s(ind_1) = exp(-k*(L + x(ind_1))) - 1 + k*(x(ind_1) + L); 
    s(ind_2) = exp(-k*(L + x(ind_2))) - 2*exp(-k*x(ind_2)) + 1 + k*(L - x(ind_2));
    s(ind_3) = exp(-k*x(ind_3)) * (1/a - 2 + a);

    s = s * oversamplingFactor / (k * k);
    
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
