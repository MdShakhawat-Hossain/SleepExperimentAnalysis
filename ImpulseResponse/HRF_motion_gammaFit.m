function [params, HRF_gamma, R2] = HRF_motion_gammaFit(HRF, Fr, method)
    % IOS_gammaFit: use gamma function to fit the HRF
    % 
    % INPUTS:
    %       HRF: numerically computed hemodynanic response function,
    %       pixcels by samples
    %       method: fitting method, 'constrained' or 'unconstrained'
    %
    % OUTPUTS:
    %       params: fitting parameters
    %       HRF_gamma: gamma function fitting
    %       R2: goodness of fit
    
    x0 = [1 1 1]; % initial guess

    % Initialize option parameters for fminsearch
    options = optimset('Display','none',...      % display no output
                       'FunValCheck','off',...  % check objective values
                       'MaxFunEvals', 1000,...  % max number of function evaluations allowed
                       'MaxIter', 1000,...      % max number of iteration allowed
                       'TolFun',1e-8,...        % termination tolerance on the function value
                       'TolX',1e-8,...          % termination tolerance on x
                       'UseParallel','always'); % always use parallel computation
    
    % ----------------------------------------- %
%     HRF = HRF(1:300); % ONLY SELECT THE FIRST 300 POINTS TO SPEEDUP PROCESS
    HRF = HRF(:);
    % ----------------------------------------- %
                   
    t = (1:length(HRF))/Fr; % time stamps for HRF
    
    if strcmp(method,'unconstrained')
        disp('IOS_gammaFit -> unconstrained fit using fminsearch');
        tic
        % unconstrained fitting
        params = fminsearch(@(x) impest_gamma(HRF,x,length(t),t), x0,options);
        toc
    end
    
    if strcmp(method, 'constrained')
        disp('IOS_gammaFit -> constrained fit using fminsearchbnd');
        tic
        % constrained fitting
        params(pix,:) = fminsearchbnd(@(x) impest_gamma(HRF,x,length(t),t),...
                                        x0, ...
                                        [-inf, realmin, realmin], ... % Lower Bound
                                        [inf,inf,inf],... % Upper Bound
                                        options);
        % W (full width at half maximum) and T (time to peak) should all be positive numbers
        toc
    end
                        
    % compute HRF estimation using gamma kernel parameters
    A = params(:,1);
    T = params(:,2);
    W = params(:,3);
    
    alpha = (T./W).^2*8.0*log(2.0);
    beta = W.^2./(T*8.0*log(2.0));

    HRF_gamma = real(A*((t/T).^alpha).*exp((t-T)/(-beta)));
    HRF_gamma = HRF_gamma(:);
    % goodness of fit
    R2 = 1-sum((HRF-HRF_gamma).^2,1)./sum(HRF.^2,1);

end

function err=impest_gamma(HRF,x,L,t)
    % impest_gamma: fit HRF using gamma kernel
    %
    % Gamma distribution kernel:
    % f(t,T,W,A) = A*(t/T).^alpha*exp((t-T)/(-beta))
    % where,
    % alpha = (T/W).^2*8.0*log(2.0)
    % beta = W^2(T*8.0*log(2.0))
    % so,
    % f(t,T,W,A) = A*(t/T).^((T/W).^2*8.0*log(2.0))*exp((t-T)/(-W^2(T*8.0*log(2.0))))
    %
    % see the following paper for more information:
    % Mariana M B Cardoso et al., The neuroimaging signal is a linear sum of
    % neurally distinct stimulus- and task-related components. Nature
    % Neuroscience, 2012.

    % Use one gamma function to fit
    A = x(1);
    T = x(2);
    W = x(3);
    
    alpha = (T/W)^2*8.0*log(2.0);
    beta = W^2/(T*8.0*log(2.0));

    HRF_gamma = A*(t/T).^alpha.*exp((t-T)/(-beta)); 
    HRF_gamma = HRF_gamma(:);

    err=sum((HRF_gamma(1:L)-HRF(1:L)).^2);
end
