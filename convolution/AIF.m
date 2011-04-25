%% AIF -- Arterial Input Function
function Cp = AIF(t)
    Cp = breastCp(t);
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
