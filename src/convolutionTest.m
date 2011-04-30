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

%     delayDemo    
    multiConvolution;
%     timingDemo
    
%     convolutionOuterLoop;
    
    '';
end

function delayDemo
    
%     y = linspace(0,0.25);
%     s = zeros(length(y));
%     for i = 1:length(y)
%         s(i) = breastCp(0.25,y(i));
%     end
%     figure, plot(y,s)
%     '';

% 
%     load fittedImages_November-04-2007_054554_2hr_20min
%     KTrans = beta(:,:,1);
%     k_ep = beta(:,:,2);
%     dt_i = 1;
%     dt_j = 1;
%     Ti = 0:28;
%     Tj = 0:28;
%     Cpi = breastCp(linspace(0,5,length(Tj)));
%     oversample_i = 1;
%     
%     
%     KTrans  = single(KTrans);
%     k_ep    = single(k_ep);
%     dt_i    = single(dt_i);
%     Ti      = int32(Ti);
%     dt_j    = single(dt_j);
%     Tj      = int32(Tj); 
%     Cpi     = single(Cpi);
%     oversample_i = single(oversample_i);
%     signal_2 = dce_mri_mex(KTrans, k_ep, dt_i, Ti, dt_j, Tj, Cpi, oversample_i);
% 
%     %figure, for t=1:29; imagesc(signal_2(:,:,t)'), axis image; drawnow; end
%     
    '';
    
    k_ep = 100;
    KTrans = 10/k_ep;

    overFactor = 5;  % Needs to be an integer
    
    oversample_i = 32;
    oversample_j = 2;
    oversample_i_b = 1 * oversample_i;
    oversample_j_b = overFactor * oversample_j;

    t0 = 0; tf = 5; T  = 50;
    %dt_j = (tf - t0) / (oversample_j*T);

    Ti = oversample_i*T;
    Tj = oversample_j*T;
    Tj2 = oversample_j_b*T;


    dt_i = (tf - t0) / Ti;
    dt_j = (tf - t0) / Tj;
    dt_j2 = (tf - t0) / Tj2;

    ti = (0:Ti-1) * dt_i;
    tj = (0:Tj-1) * dt_j;
    tj2 = (0:Tj2-1) * dt_j2;

    tj_vec = single(0:Tj-1) * dt_j;
    tj_vec = permute(tj_vec, [1 3 2]);  % Make it 1x1xTj
    
    
    % Init Cpi vector
    Cpi = breastCp(ti);
    [v,ind] = max(Cpi)
    f = @(x) 5 - breastCp(x);
    ti_min = fminsearch(f, ti(ind));
    [v,ind] = min(abs(ti - ti_min));
    [ti(ind), Cpi(ind)]
    
    
    
    h2 = figure;
    hold all
    title('Difference')
    
    h1 = figure;
    hold all
    %plot(ti, Cpi*10, 'LineWidth', 3, 'Color', [.8 .8 1])
    
    
%     Ti = oversample_i*T;
%     Tj = oversample_j*T;
%     
%     dt_i = (tf - t0) / Ti;
%     dt_j = (tf - t0) / Tj;
    
    delays = linspace(0, dt_j, 10);
    for k = 1:length(delays)
        delay = delays(k);

        signal_a = convolutionOversampleRun (KTrans, k_ep, delay, oversample_i, oversample_j);
        signal_b = convolutionOversampleRun (KTrans, k_ep, delay, oversample_i_b, oversample_j_b);

        
        v = (k-1) / max(1, length(delays)-1);
        color_a = [1-v,   v, v];
        color_b = [0.9*(1-v), 0.75*(1-v), v];

        s_a = real(signal_a);
        s_a = permute(s_a, [3 1 2]);
        s_a = s_a(:,:);
        t_a = tj(:);

        s_b = real(signal_b);
        s_b = permute(s_b, [3 1 2]);
        s_b = s_b(:,:);
        t_b = tj2(:);
        
        s_c = s_b(1:overFactor:end);
        s_c = s_c(1:min(end,length(s_a)));  % Force to be same size if something went wrong
        t_c = t_a;
        
        set(0,'CurrentFigure', h1)
        plot(t_a, s_a, '--', 'Color', color_a) 
        plot(t_b, s_b,  '-', 'Color', color_b) 
        plot(t_c, s_c,  'o', 'Color', color_b)

        set(0,'CurrentFigure', h2)
        plot(t_a, s_a - s_c, 'Color', color_a)
        
    end
    
    '';
end


function timingDemo
% 
%     load fittedImages_November-04-2007_054554_2hr_20min
%     KTrans = beta(:,:,1);
%     k_ep = beta(:,:,2);
%     dt_i = 1;
%     dt_j = 1;
%     Ti = 0:28;
%     Tj = 0:28;
%     Cpi = breastCp(linspace(0,5,length(Tj)));
%     oversample_i = 1;
%     
%     
%     KTrans  = single(KTrans);
%     k_ep    = single(k_ep);
%     dt_i    = single(dt_i);
%     Ti      = int32(Ti);
%     dt_j    = single(dt_j);
%     Tj      = int32(Tj); 
%     Cpi     = single(Cpi);
%     oversample_i = single(oversample_i);
%     signal_2 = dce_mri_mex(KTrans, k_ep, dt_i, Ti, dt_j, Tj, Cpi, oversample_i);
% 
%     %figure, for t=1:29; imagesc(signal_2(:,:,t)'), axis image; drawnow; end
%     
    '';
    %convolutionOuterLoop;
    %multiConvolution;
    %[time_matlab, time_cuda, RMSE, nRMSE] = convolutionRun(16,16,100)
    
    Xs = 2.^(0:6);
    Ys = 2.^(0:6);
    %Ts = 2.^[5:7];
    Ts = 64;
    time_matlab = zeros(length(Xs),length(Ts));
    time_cuda   = zeros(length(Xs),length(Ts));
    RMSE        = zeros(length(Xs),length(Ts));
    nRMSE       = zeros(length(Xs),length(Ts));
%     for X = Xs
%         for Y = Ys
%             for T = Ts
    for i = 1:length(Xs)
        %for j = 1:length(Ys)
            k = 1;
            X = Xs(i);
            Y = Ys(i);
            T = Ts(k);
            %for T = Ts
                [time_matlab(i,k), time_cuda(i,k), RMSE(i,k), nRMSE(i,k)] = convolutionRun(X,Y,T);
            %end
        %end
    end

    time_matlab
    time_cuda
    RMSE
    nRMSE

    displayInUnits = @(src,evt) set(gca,...
        'XTickLabel', cellfun(@num2str, num2cell(get(gca, 'XTick')), 'UniformOutput',false), ...
        'YTickLabel', cellfun(@num2str, num2cell(get(gca, 'YTick')), 'UniformOutput',false) ...
    );
    
    %set(gcf,'ResizeFcn', displayInUnits)
    
    figure('ResizeFcn', displayInUnits)
    %hold all
    fontsize_title  = 48;
    fontsize        = floor(3/4*fontsize_title);
    fontsize_axis   = floor(1/3*fontsize_title);
    markersize = 25;
    linewidth = 5;
    loglog(Xs.*Ys, time_matlab, '--X', Xs.*Ys, time_cuda, '-*', 'LineWidth', linewidth, 'MarkerSize', markersize)
    grid minor
    title('Time in seconds as a function of problem size (X x Y)', 'FontSize', fontsize_title)
    %legend('\fontsize{32}Matlab', '\fontsize{32}Cuda')
    legend(sprintf('\\fontsize{%d}CPU',fontsize_axis), sprintf('\\fontsize{%d}GPU',fontsize_axis))
    %legend('CPU', 'GPU')
    ylabel('t (sec)', 'Fontsize', fontsize)
    xlabel('Number of Voxels', 'Fontsize', fontsize)
    set(gca,'FontSize', fontsize_axis)
    set(gca,'YTickLabel', cellfun(@num2str, num2cell(get(gca, 'YTick')), 'UniformOutput',false))
    set(gca,'XTickLabel', cellfun(@num2str, num2cell(get(gca, 'XTick')), 'UniformOutput',false))
   
    %loglog(Xs.*Ys, time_cuda, '-o')
    %diag(time_matlab)
    %diag(time_cuda)
    '';
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
    % Init time variables
    t0 = 0; 
    tf = 5;
    T  = 50;

    oversample_i = 16;
    oversample_j = 4;
    
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

    k_ep = 100;
    KTrans = k_ep;
    
    
    % Convert to 32-bit floats for everything inputted into dce_mri_mex()
    KTrans  = single(KTrans);
    k_ep    = single(k_ep);
    dt_i    = single(dt_i);
    dt_j    = single(dt_j);
    Ti      = int32(Ti);
    Tj      = int32(Tj); 
    Cpi     = single(Cpi);
    oversample_i = single(oversample_i);
    oversample_j = single(oversample_j);


    tj_vec = single(0:Tj-1) * dt_j;
    tj_vec = permute(tj_vec, [1 3 2]);  % Make it 1x1xTj
    

    h2 = figure;
    hold all
    title('Difference')
    
    h1 = figure;
    hold all
    plot(ti, Cpi*10, 'LineWidth', 3, 'Color', [.8 .8 1])
    
    %for k_ep = logspace(log10(0.01), log10(1000), 20);
    %for delay = delays
    
    delays = linspace(0, dt_j, 10);
    for k = 1:length(delays)
        delay = delays(k);

        signal = zeros(1,1,Tj,'single');
        for i = 1:Ti
            signal = signal + Cpi(i) * convolutionFromMapleMatrized(tj_vec - delay, k_ep, ti(i), oversample_i);
        end
        signal = bsxfun(@times,signal,KTrans);

        signal_alt = zeros(1,1,Tj,'single');
        for i = 1:Ti
            signal_alt = signal_alt + Cpi(i) * convolutionFromMapleMatrized(tj_vec, k_ep, ti(i) + delay, oversample_i);
        end
        signal_alt = bsxfun(@times,signal_alt,KTrans);

%         signal_alt = dce_mri_mex(single(KTrans), single(k_ep), single(dt_i), int32(Ti), single(dt_j), int32(Tj), single(Cpi), single(oversample_i));
        
        v = (k-1) / max(1, length(delays)-1);
        color = [v, 1-v, 1-v];
        
        set(0,'CurrentFigure', h1)
        plot(tj, real(signal(:)), 'Color', color)

        set(0,'CurrentFigure', h2)
        plot(tj, real(signal(:) - signal_alt(:)), 'Color', color)
        
    end
    '';
end

%%
function signal = convolutionOversampleRun (KTrans, k_ep, delay, oversample_i, oversample_j)
    [X,Y] = size(KTrans);
    
    % Init time variables
    t0 = 0; 
    tf = 5;
    T  = 50;

%     oversample_i = 16;
%     oversample_j = 4;
    
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

%     k_ep = 100;
%     KTrans = k_ep;
    
    
    % Convert to 32-bit floats for everything inputted into dce_mri_mex()
    KTrans  = single(KTrans);
    k_ep    = single(k_ep);
    dt_i    = single(dt_i);
    dt_j    = single(dt_j);
    Ti      = int32(Ti);
    Tj      = int32(Tj); 
    Cpi     = single(Cpi);
    oversample_i = single(oversample_i);
    oversample_j = single(oversample_j);


    tj_vec = single(0:Tj-1) * dt_j;
    tj_vec = permute(tj_vec, [1 3 2]);  % Make it 1x1xTj
    

%     h2 = figure;
%     hold all
%     title('Difference')
%     
%     h1 = figure;
%     hold all
%     plot(ti, Cpi*10, 'LineWidth', 3, 'Color', [.8 .8 1])
    
    %for k_ep = logspace(log10(0.01), log10(1000), 20);
    %for delay = delays
    
    %delays = linspace(0, dt_j, 10);
%     for k = 1:length(delays)
%         delay = delays(k);
%         
%         sj = zeros(1,Tj);
%         for i=1:Ti
%             signal_part = convolutionInnerLoop(tj, k_ep, ti(i), oversample_i);
%             sj = sj + Cpi(i) * KTrans * signal_part;
%         end

        %sj = zeros(X,Y,Tj,'single');
        signal = zeros(X,Y,Tj,'single');
        for i = 1:Ti
            signal = signal + Cpi(i) * convolutionFromMapleMatrized(tj_vec - delay, k_ep, ti(i), oversample_i);
        end
        signal = bsxfun(@times,signal,KTrans);

%         signal_alt = zeros(1,1,Tj,'single');
%         for i = 1:Ti
%             signal_alt = signal_alt + Cpi(i) * convolutionFromMapleMatrized(tj_vec, k_ep, ti(i) + delay, oversample_i);
%         end
%         signal_alt = bsxfun(@times,signal_alt,KTrans);

%         signal_alt = dce_mri_mex(single(KTrans), single(k_ep), single(dt_i), int32(Ti), single(dt_j), int32(Tj), single(Cpi), single(oversample_i));
        
%         v = (k-1) / max(1, length(delays)-1);
%         color = [v, 1-v, 1-v];
        
%         set(0,'CurrentFigure', h1)
%         plot(tj, real(signal(:)), 'Color', color)

%         set(0,'CurrentFigure', h2)
%         plot(tj, real(signal(:) - signal_alt(:)), 'Color', color)
        
%     end
    '';
end


%%
function [time_matlab, time_cuda, RMSE, nRMSE] = convolutionRun(X,Y,T)
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
    
    [mX,mY] = meshgrid(linspace(-3,3,X),linspace(-3,3,Y)); Z = peaks(mX,mY)'; 
    b = 1; c = 4; f = @(x) 1/b * x./(1 + (x*c-c).^2); %plot(x, f(x))
    Z = f(Z);
    Z = Z * 10;
    Z = max(Z, 1e-1);
        
    k_ep = Z;
    KTrans = Z/10;

    % Init time variables
    t0 = 0; 
    tf = 5;

    oversample_i = 8;
    oversample_j = 2;
    
    Ti = oversample_i*T;
    Tj = oversample_j*T;
    
    dt_i = (tf - t0) / Ti;
    dt_j = (tf - t0) / Tj;
    

    ti = (0:Ti-1) * dt_i;
    tj = (0:Tj-1) * dt_j;


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
    

    fprintf('Beginning the timings\n')

    
    % Reference implementation
	tj_vec = single(0:Tj-1) * dt_j; 
    tj_vec = permute(tj_vec, [1 3 2]);  % Make it 1x1xTj
    
    time_matlab = tic;
    sj = zeros(X,Y,Tj,'single');
    for i = 1:Ti
        sj = sj + Cpi(i) * convolutionFromMapleMatrized(tj_vec, k_ep, ti(i), oversample_i);
    end
    signal_0 = bsxfun(@times,sj,KTrans);
    time_matlab = toc(time_matlab);

    
    
    
    
    % Matrix version of mex-wrapper call
    time_cuda = tic;
    signal_2 = dce_mri_mex(KTrans, k_ep, dt_i, Ti, dt_j, Tj, Cpi, oversample_i);
    time_cuda = toc(time_cuda);
    
%     signal_1 = real(signal_1);
    signal_2 = real(signal_2);

    

    % Compare cuda versus matlab
    RMSE  = norm(signal_2(:) - signal_0(:)) / sqrt(X*Y*double(Tj));
    nRMSE = RMSE / norm(signal_0(:));
    
%     fprintf('Cuda versus matlab:  RMSE: %g\n',  RMSE)
%     fprintf('Cuda versus matlab: nRMSE: %g\n', nRMSE)
        
end






%%
function signal = multiConvolution()
    % Image dimensions
%     X = 64; Y = 64;
    X = 16; Y = 16;
%     X = 8;  Y = 8;

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
    %Z = max(Z, 1e-1);
    Z = max(Z, 1/8);
    
    %figure, imagesc(Z'), axis image; colorbar
    
    k_ep = Z;
    KTrans = Z/10;

    
    figure, imagesc(Z), axis image; colorbar, title('k_ep')
    
    
    % Init time variables
    t0 = 0; 
    tf = 5;
    T  = 50;

    oversample_i = 4;
    oversample_j = 16;
    
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
%     signal_0 = zeros(X,Y,Tj,'single');
	tj_vec = single(0:Tj-1) * dt_j; 
    tj_vec = permute(tj_vec, [1 3 2]);  % Make it 1x1xTj
    tic
%     for x = 1:X
%         for y = 1:Y
%             sj = zeros(1,1,Tj,'single');
%             for i = 1:Ti
%                 %si = si + Cpj(j) * KTrans * convolutionFromMapleVectorized(ti, k_ep, tj(j), oversample_j);
%                 sj = sj + Cpi(i) * KTrans(x,y) * convolutionFromMapleVectorized(tj_vec, k_ep(x,y), ti(i), oversample_i);
%                 %signal_0(x,y,:) = dce_mri_mex( KTrans(x,y), k_ep(x,y), dt_i, Ti, dt_j, Tj, Cpi, oversample_i );
%             end
%             signal_0(x,y,:) = sj;
%         end
%     end
    sj = zeros(X,Y,Tj,'single');
    for i = 1:Ti
        %foo = convolutionFromMapleMatrized(tj_vec, k_ep, ti(i), oversample_i);
        sj = sj + Cpi(i) * convolutionFromMapleMatrized(tj_vec, k_ep, ti(i), oversample_i);
        %s = convolutionFromMapleMatrized(t, k, t_0, oversamplingFactor)
    end
    signal_0 = bsxfun(@times,sj,KTrans);
    fprintf('Reference implementation (host): '), toc
   % norm(signal_0(:) - signal_4(:))
    '';
    
    
    
%     % Scalar version of mex-wrapper call
%     signal_1 = zeros(X,Y,Tj,'single');
%     tic
%     for x = 1:X
%         for y = 1:Y
%             signal_1(x,y,:) = dce_mri_mex( KTrans(x,y), k_ep(x,y), dt_i, Ti, dt_j, Tj, Cpi, oversample_i );
%         end
%     end
%     fprintf('Scalar version of mex-wrapper call: '), toc


    d = 2^10;
    KTrans = round(d*KTrans)/d;
    k_ep = round(d*k_ep)/d;
    Cpi = round(d*Cpi)/d;
    

    % Matrix version of mex-wrapper call
    tic
    signal_3 = dce_mri_mex(KTrans, k_ep, dt_i, Ti, dt_j, Tj, Cpi, oversample_i, false);
    fprintf('Matrix version of mex-wrapper call (CPU): '), toc

    
    % Matrix version of mex-wrapper call
    tic
    signal_2 = dce_mri_mex(KTrans, k_ep, dt_i, Ti, dt_j, Tj, Cpi, oversample_i, true);
    fprintf('Matrix version of mex-wrapper call (GPU): '), toc



    %colors = rand(X*Y,3);
%     [v,ind] = sort(k_ep(:));
%     colors(ind,1:3) = jet(X*Y);
    
    u = log10(k_ep(:));
    u = (u - min(u)) / (max(u) - min(u) + eps);
    
    palette = jet(1000);
    ind = min(size(palette,1), max(1, 1 + floor(u*max(1, size(palette,1) - 1))));
    colors(:,1:3) = palette(ind,1:3);
    
    
    %cdef = 'white';
    cdef = 'black';
    figure, colordef(gcf, cdef), set(gca, 'ColorOrder', colors), hold all; plot(tj, reshape(permute(real(signal_2),[3 1 2]), [Tj X*Y]))
    figure, colordef(gcf, cdef), set(gca, 'ColorOrder', colors), hold all; plot(tj, reshape(permute(real(signal_3),[3 1 2]), [Tj X*Y]))
    figure, colordef(gcf, cdef), set(gca, 'ColorOrder', colors), hold all; plot(tj, reshape(permute(real(signal_2),[3 1 2]), [Tj X*Y]), tj, reshape(permute(real(signal_3),[3 1 2]), [Tj X*Y]), '--')
    figure, colordef(gcf, cdef), set(gca, 'ColorOrder', colors), hold all; plot(tj, reshape(permute(real(signal_2),[3 1 2]), [Tj X*Y]) - reshape(permute(real(signal_3),[3 1 2]), [Tj X*Y]))
    
    figure, colordef(gcf, cdef), set(gca, 'ColorOrder', colors), hold all; plot(tj, reshape(permute(real(signal_2),[3 1 2]), [Tj X*Y]) - reshape(permute(real(signal_0),[3 1 2]), [Tj X*Y]))
    figure, colordef(gcf, cdef), set(gca, 'ColorOrder', colors), hold all; plot(tj, reshape(permute(real(signal_3),[3 1 2]), [Tj X*Y]) - reshape(permute(real(signal_0),[3 1 2]), [Tj X*Y]))
    %figure, colordef(gcf, cdef), set(gca, 'ColorOrder', colors), hold all; plot(tj, reshape(permute(real(signal_2),[3 1 2]), [Tj X*Y]), tj, reshape(permute(real(signal_3),[3 1 2]), [Tj X*Y]), '--', tj, reshape(permute(real(signal_0),[3 1 2]), [Tj X*Y]), ':')
    
    
%     signal_1 = real(signal_1);
    signal_2 = real(signal_2);
    
    signal = signal_2;  % Final return value of this function

%     
%     % Compare matrix versus scalar calls
%     RMSE  = norm(signal_1(:) - signal_2(:)) / sqrt(X*Y*double(Tj));
%     nRMSE = RMSE / norm(signal_1(:));
%     
%     fprintf('Matrix versus scalar:  RMSE: %g\n',  RMSE)
%     fprintf('Matrix versus scalar: nRMSE: %g\n', nRMSE)
    

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
function s = convolutionFromMapleMatrized(t, k, t_0, oversamplingFactor)
    x = t - t_0;
    L = 1/oversamplingFactor;
    s = zeros([size(k),length(t)]);

    
    ind_1 = (x > -L  &  x <= 0);
    ind_2 = (x >  0  &  x <= L);
    ind_3 = (x > L);
    
    kx1 = bsxfun(@times, k, x(ind_1));
    kx2 = bsxfun(@times, k, x(ind_2));
    kx3 = bsxfun(@times, k, x(ind_3));
    
    a = exp(k*L);
    %s(ind_1) = exp(-k.*(L + x(ind_1))) - 1 + k.*(x(ind_1) + L); 
    s(:,:,ind_1) = bsxfun(@plus, exp(-bsxfun(@times, k, L + x(ind_1))) - 1 + kx1, k*L); 
    %s(ind_2) = exp(-k.*(L + x(ind_2))) - 2*exp(-k.*x(ind_2)) + 1 + k.*(L - x(ind_2));
    s(:,:,ind_2) = bsxfun(@plus, exp(-bsxfun(@times, k, (L + x(ind_2)))) - 2*exp(-kx2) + 1 - kx2,  k*L);
    s(:,:,ind_3) = bsxfun(@times, exp(-kx3), (1./a - 2 + a));

    s = bsxfun(@rdivide, oversamplingFactor*s, (k .* k));
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
    if nargin < 2, dt = 0.25; end;
    %if (nargin < 2) dt = 32/60; end;

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
