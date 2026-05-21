function [A_opt, out] = NewtonSolverReparameterization(lambda_true, weights, F, dE, A_init, opts)
% minimize_KL_Newton: Minimize D_JS(Pm, Q(A)) where Q(A) has Poisson rates
% tilde_lambda_i(A) = ∫ w_i(E) * exp(-sum_j A_j f_j(E)) dE  (discretized).
%
%   Inputs:
%   lambda_true : [M x 1] ground-truth rates (fixed)
%   weights     : [M x nE] w_i(E)
%   F           : [N x nE] basis f_j(E)
%   dE          : scalar or [1 x nE] energy widths
%   A_init      : [N x 1] initial guess for A (nonnegative)
%   t           : scalar t-value in log-barrier optimization
%   opts        : struct with fields (all optional)
%       maxIter (default 100)
%       gtol    (default 1e-6)   gradient norm tolerance
%       xtol    (default 1e-8)   parameter step tolerance
%       ftol    (default 1e-10)  objective decrease tolerance
%       mu0     (default 1e-2)   initial LM damping
%       floorRate (default 1e-12)
%       verbose (default 1)
%
% Outputs:
%   A_opt : [N x 1] optimized A
%   out   : struct with diagnostics (history of f, grad norms, etc.)

    if nargin < 6, opts = struct(); end
    getopt = @(name,def) (isfield(opts,name) * opts.(name) + ~isfield(opts,name) * def);

    maxIter   = getopt('maxIter',   100);
    gtol      = getopt('gtol',      1e-6);
    xtol      = getopt('xtol',      1e-8);
    ftol      = getopt('ftol',      1e-10);
    mu        = getopt('mu0',       1e-3);
    floorRate = getopt('floorRate', 1e-12);
    verbose   = getopt('verbose',   1);
    mu_reg    = getopt('mu_reg',     1e-8);

    % Ensure shapes
    [M, nE] = size(weights);
    [N, nE2]= size(F);
    if nE ~= nE2, error('weights and F must have same # of energy bins'); end
    if isscalar(dE), dE = dE * ones(1, nE); end
    if ~iscolumn(lambda_true), lambda_true = lambda_true(:); end
    if ~iscolumn(A_init), A_init = A_init(:); end

    % Reparameterize: A = exp(B)  (enforce A >= 0)
    A_init = max(A_init, floorRate);

    hist.f = zeros(maxIter+1,1);
    hist.gnorm = zeros(maxIter+1,1);
    
    %% Initial evaluation
    
    B = log(max(A_init,1e-12));
    A = exp(B);
    
    [~, gA, HA, ~, ~, ~] = NewtonStep(A, lambda_true, weights, F, dE, floorRate);
    
    %% Scaling damping according to the problem. 
    
    HB0 = diag(A) * HA * diag(A) + diag(A .* gA);
    HB0 = 0.5*(HB0 + HB0.');
    Hscale = max(1, max(abs(diag(HB0))));
    mu = opts.mu0 * Hscale;       

    for k=1:maxIter
        [f, gA, HA, ~, ~, ~] = NewtonStep(A, lambda_true, weights, F, dE, floorRate);
        
        %% Reparameterization A = exp(B).
        
        gB = A .* gA;  % Gradient
        HB = diag(A) * HA * diag(A) + diag(gA .* A); % Hessian
        
        %% regularize Hessian
        
        HB = 0.5*(HB + HB.'); % maintain symmetry. Safeguards for numerical errors.
        
        
        Hscale = max(1, max(abs(diag(HB))));
        reg = max(mu + mu_reg, 1e-12 * Hscale);
        Hreg = HB + reg * eye(N);
        
        %% Compute Newton direction:
        % We try Cholesky-based factorization to ensure numerical stability when Hessian has small norm
        
        attempts = 0;
        jitter = reg;
        while attempts < 10
            try
                R = chol(Hreg);                  % will error if not SPD
                pB = - R \ (R' \ gB);            % efficient SPD solve
                break;
            catch
                % add more diagonal jitter until we succeed
                jitter = max(1e-6, 10 * jitter);
                Hreg = HB + jitter * eye(N);
                attempts = attempts + 1;
            end
        end
        if attempts == 10
            warning('Could not make Hessian SPD after jitter attempts; using backslash as fallback.');
            pB = - Hreg \ gB;  % fallback
        end
                

        %% stopping on gradient norm
        
        if norm(gB) < gtol
            if verbose, fprintf('Converged: grad norm %.3e < gtol\n', norm(gB)); end
            break;
        end
        
        %% Armijo backtracking in B-space
        
        c1 = 1e-4; % Armijo paramter
        [alpha, f_trial] = ArmijoBacktracking(); % Armijo etc using f(B)->f(A=exp(B))
        %% Update damping with Levenberg-Marquardt rule
        
        LevenbergMarquardtUpdate();
        
        %% Accept step
        
        B_new = B + alpha * pB;
        step_norm = norm(alpha * pB);
        B = B_new;
        f = f_trial;           % update current objective to accepted trial
        A = exp(B);
        %% store history

        hist_f(k+1) = f;
        hist_gnorm(k+1) = norm(gB);
        hist_step(k+1) = norm(step_norm);   
        hist_A = A;

        if verbose
            fprintf('Iter %3d: f=% .6e, ||gB||=%.3e, step=%.3e, alpha=%.3e \n', k, f, hist_gnorm(k+1), step_norm, alpha);
        end

        %% stopping on step size
        if step_norm < xtol * (1 + norm(B))
            if verbose, fprintf('Converged: step small\n'); end
            break;
        end
        

        
        
    end
   A_opt = exp(B);
    % Trim history arrays
    last = min(k+1, length(hist_f));
    out.hist.fval  = hist_f(1:last);
    out.hist.gnorm = hist_gnorm(1:last);
    out.hist.step  = hist_step(1:last);
    out.iter = k;
    out.final_f = f;
    out.final_grad_norm = norm(gB);
    out.A_opt = A_opt;


function [alpha, f_trial] = ArmijoBacktracking()
    
    alpha = 1;
    f_curr = f;
    B_curr = B;
    gB_curr = gB;
    maxbacktracks = 20;
   
    for i=1:maxbacktracks
        B_trial = B_curr + alpha * pB;
        A_trial = exp(B_trial);
        
        [f_trial,~, ~, ~, ~, ~] = NewtonStep(A_trial, lambda_true, weights, F, dE, floorRate);
        
        % Armijo condition (gA is gradient wrt A)
        if f_trial <= f_curr + c1 * alpha * (gB_curr.' * pB)
            % sufficient decrease by Armijo condition
            return;
        end
        
        % reduce step
        alpha = 0.5 * alpha;
       

    end
    
    % if we exit the loop without acceptance, return smallest alpha
    if alpha < 1e-12
        alpha = 0;
    end
end

function LevenbergMarquardtUpdate()
    
    %% LM -params
    
    rho_up = 0.75;    % threshold for good step
    rho_down = 0.25;  % threshold for poor step
    factor_up = 0.1;  % shrink mu by this factor if good
    factor_down = 10; % increase mu by this factor if poor
    mu_min = 1e-16;
    mu_max = 1e16;
    
    %% Compute step quality
    
    s = pB * alpha; % Acutal step
    
    actual_reduction = f - f_trial;
    predicted_reduction = - gB' * s - 0.5 * s' * HB * s;
    
    eps = 1e-10;
    predicted_reduction = max(eps, predicted_reduction);
    
    rho = actual_reduction / predicted_reduction;
    
    %% Update mu based on rho
    
    if rho > rho_up
        mu = max(mu * factor_up, mu_min);
    elseif rho < rho_down
        mu = min(mu * factor_down, mu_max);
    end
    
end

end

function [f, gA, HA, lambda_model, J, A] = NewtonStep(A, lambda_true, weights, F, dE, floorRate)
    % Compute value, gradient (w.r.t A), and Hessian (w.r.t A) of
    
    useKL = 0; 
    useKLsym = 1;
    
    % Shapes
    [M, nE] = size(weights);
    [N, ~]  = size(F);
    A = A(:);

    % Exponent term: z(E) = A' * F(:,E), then exp(-z)
    z = A.' * F;                  % [1 x nE]
    exp_term = exp(-z);           % [1 x nE]
    % Discretized integrals
    % S = weights .* (replicate exp_term) .* dE
    S = weights .* (ones(M,1) * (exp_term .* dE));   % [M x nE]

    % Rates
    lambda_model = sum(S, 2);                         % [M x 1]
    lambda_model = max(lambda_model, floorRate);      % positivity floor
    
    % Jacobian of tilde_lambda w.r.t A: J = - S * F'   (M x N)
    J = - (S * F.');                                  % [M x N]
    
    if ~useKLsym % Minimizes KL-divergence: KL(P || Q)
        % Residual weight vector
        r = 1 - (lambda_true ./ lambda_model);           % [M x 1]

        % Value
        f = sum( lambda_true .* (log(lambda_true ./ lambda_model)) - (lambda_true - lambda_model));
        f = abs(f);
        
        % Hessian w.r.t A:
        % GN term: J' * diag(w) * J with w_i = lambda_i / tilde_lambda_i^2
        wgn = lambda_true ./ (lambda_model.^2);           % [M x 1]

    else % Minimizes symmetrised KL-divergece: KL(P || Q) + KL(Q || P).
        % Residual weight vector
        r = 1 - (lambda_true ./ lambda_model) + log(lambda_model ./ lambda_true);           % [M x 1]

        % Value
        f = sum( lambda_true .* (log(lambda_true ./ lambda_model)) + lambda_model .* (log(lambda_model ./ lambda_true)));
        f = abs(f); % Ensure f is positive for numerical stability.
        
        
        % Hessian w.r.t A:
        % GN term: J' * diag(w) * J with w_i = lambda_i / tilde_lambda_i^2
        wgn = (lambda_true + lambda_model) ./ (lambda_model.^2);           % [M x 1]
    end
    
    % Gradient w.r.t A: gA = J' * r
    gA = J.' * r;                                     % [N x 1]
    
    Hgn = J.' * (bsxfun(@times, J, wgn));            % [N x N]
    % Curvature term: sum_i (1 - lambda_i/tilde_lambda_i) * Htilde_i
    % Htilde_i = ∫ w_i(E) f(E) f(E)^T exp(-A·f(E)) dE
    % Efficient sum: F * diag(t) * F', where t(E) = (c' * w(:,E)) * exp_term(E) * dE(E)
    c = r.' * weights;                                % [1 x nE]
    Hcurv = (F .* c .* exp_term .* dE) * F.'; % [N x N]; F*diag(t)*F' == (F.*t)*F'

    HA = (Hgn + Hcurv);
    
    % Numerical guarding (optional, keeps symmetry)
    
    HA = 0.5*(HA + HA.');
end