function [A_opt, out] = NewtonSolver(lambda_true, weights, F, dE, A_init, opts)
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
    A = max(A_init, floorRate);

    hist.f = zeros(maxIter+1,1);
    hist.gnorm = zeros(maxIter+1,1);
    
    %% Initial evaluation
    

    
    [~, gA, HA, ~, ~, ~] = NewtonStep(A, lambda_true, weights, F, dE, floorRate);
    
    %% Scaling damping according to the problem. 

    Hscale = max(1, max(abs(diag(HA))));
    mu = opts.mu0 * Hscale;       

    for k=1:maxIter
        [f, gA, HA, ~, ~, ~] = NewtonStep(A, lambda_true, weights, F, dE, floorRate);
        
        
        %% regularize Hessian
        
        HA = 0.5*(HA + HA.'); % maintain symmetry. Safeguards for numerical errors.
        
        
        Hscale = max(1, max(abs(diag(HA))));
        reg = max(mu + mu_reg, 1e-12 * Hscale);
        Hreg = HA + reg * eye(N);
        
        %% Compute Newton direction:
        % We try Cholesky-based factorization to ensure numerical stability when Hessian has small norm
        
        attempts = 0;
        jitter = reg;
        while attempts < 10
            try
                R = chol(Hreg);                  % will error if not SPD
                p = - R \ (R' \ gA);            % efficient SPD solve
                break;
            catch
                % add more diagonal jitter until we succeed
                jitter = max(1e-6, 10 * jitter);
                Hreg = HA + jitter * eye(N);
                attempts = attempts + 1;
            end
        end
        if attempts == 10
            warning('Could not make Hessian SPD after jitter attempts; using backslash as fallback.');
            p = - Hreg \ gA;  % fallback
        end
                

        %% stopping on gradient norm
        
        if norm(gA) < gtol
            if verbose, fprintf('Converged: grad norm %.3e < gtol\n', norm(gA)); end
            break;
        end
        
        %% Armijo backtracking in B-space
        
        c1 = 1e-4; % Armijo paramter
        [alpha, f_trial] = ArmijoBacktracking(); % Armijo etc using f(B)->f(A=exp(B))
        
        %% Update damping with Levenberg-Marquardt rule
        
        LevenbergMarquardtUpdate();
        
        %% Accept step
        
        A_new = A + alpha * p;
        step_norm = norm(alpha * p);
        A = A_new;
        f = f_trial;           % update current objective to accepted trial

        %% store history

        hist_f(k+1) = f;
        hist_gnorm(k+1) = norm(gA);
        hist_step(k+1) = norm(step_norm);   
        hist_A = A;

        if verbose
            fprintf('Iter %3d: f=% .6e, ||gB||=%.3e, step=%.3e, alpha=%.3e \n', k, f, hist_gnorm(k+1), step_norm, alpha);
        end

        %% stopping on step size
        if step_norm < xtol * (1 + norm(A))
            if verbose, fprintf('Converged: step small\n'); end
            break;
        end
        

        
        
    end
    A_opt = A;
    % Trim history arrays
    last = min(k+1, length(hist_f));
    out.hist.fval  = hist_f(1:last);
    out.hist.gnorm = hist_gnorm(1:last);
    out.hist.step  = hist_step(1:last);
    out.iter = k;
    out.final_f = f;
    out.final_grad_norm = norm(gA);
    out.A_opt = A_opt;


function [alpha, f_trial] = ArmijoBacktracking()
    
    alpha = 1;
    f_curr = f;
    A_curr = A;
    gA_curr = gA;
    maxbacktracks = 20;
   
    for i=1:maxbacktracks
        A_trial = A_curr + alpha * p;
        
        [f_trial,~, ~, ~, ~, ~] = NewtonStep(A_trial, lambda_true, weights, F, dE, floorRate);
        
        % Armijo condition (gA is gradient wrt A)
        if f_trial <= f_curr + c1 * alpha * (gA_curr.' * p)
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
    
    s = p * alpha; % Acutal step
    
    actual_reduction = f - f_trial;
    predicted_reduction = - gA' * s - 0.5 * s' * HA * s;
    
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




function [f, gA, HA] = NewtonStepBatch(A, lambda_true, weights, F, dE, nMat, nBins, nBatch, floorRate)
    % Robust, memory-safe batch version consistent with NewtonStep (scalar).
    % Input shapes expected:
    %   A          : 1 x nMat x 1 x 1 x nBatch   (as you build it)
    %   lambda_true: 1 x 1 x nBins x 1 x nBatch  (as produced earlier)
    %   weights    : either 1 x 1 x nBins x nE  OR  nBins x nE
    %   F          : nMat x nE
    %   dE         : scalar or [1 x nE]
    %
    % Outputs:
    %   f  : [nBatch x 1]
    %   gA : [nMat x nBatch]
    %   HA : [nMat x nMat x nBatch]

    % ---- unify weights to 2-D (nBins x nE) ----
    
    if ndims(weights) == 4
        w = squeeze(weights(1,1,:,:));   % -> [nBins x nE]
    elseif ismatrix(weights)
        w = weights;
    else
        error('Unexpected weights shape in NewtonStepBatch.');
    end
    [wBins, wE] = size(w);
    if wBins ~= nBins || wE ~= size(F,2)
        error('weights size does not match nBins/nE.');
    end

    % ---- unify lambda_true to [nBins x nBatch] ----
    lambda_true2 = squeeze(lambda_true);  % should become [nBins x nBatch]
    if isvector(lambda_true2)
        lambda_true2 = reshape(lambda_true2, [nBins, nBatch]);
    end
    if any(size(lambda_true2) ~= [nBins, nBatch])
        error('lambda_true has unexpected shape after squeeze.');
    end

    % ---- A -> [nMat x nBatch] and forward projection ----
    A_mat = reshape(A, [nMat, nBatch]);    % [nMat x nBatch]
    Z = A_mat.' * F;                       % [nBatch x nE]
    exp_term = exp(-Z);                    % [nBatch x nE]

    % ---- ensure dE is 1 x nE ----
    if isscalar(dE)
        dE_row = dE;
    else
        dE_row = reshape(dE, [1, numel(dE)]);
    end

    % ---- predicted lambdas: lambda_contrast (nBins x nBatch) ----
    % X: (nBatch x nE) .* (1 x nE) -> nBatch x nE, then multiply w * X'
    X = bsxfun(@times, exp_term, dE_row);      % [nBatch x nE]
    lambda_contrast = w * X.';                 % (nBins x nE) * (nE x nBatch) -> nBins x nBatch
    lambda_contrast = max(lambda_contrast, floorRate);

    % ---- residual r and objective f ----
    r = 1 - lambda_true2 ./ lambda_contrast + log(lambda_contrast ./ lambda_true2);   % [nBins x nBatch]
    fvec = sum(lambda_true2 .* log(lambda_true2 ./ lambda_contrast) + ...
               lambda_contrast .* log(lambda_contrast ./ lambda_true2), 1);          % [1 x nBatch]
    f = fvec(:);  % [nBatch x 1]

    % ---- prepare outputs ----
    gA = zeros(nMat, nBatch);
    HA = zeros(nMat, nMat, nBatch);

    % ---- loop over batches (memory-safe) ----
    for b = 1:nBatch
        expb = exp_term(b,:);                     % [1 x nE]
        tb = bsxfun(@times, expb, dE_row);        % [1 x nE]  (exp* dE)

        % S_b = w .* tb  (broadcast tb across rows)
        Sb = bsxfun(@times, w, tb);               % [nBins x nE]

        % Jacobian Jb = - S_b * F.'  -> [nBins x nMat]
        Jb = - (Sb * F.');                        % [nBins x nMat]

        rb = r(:, b);                             % [nBins x 1]

        % gradient: gA(:,b) = Jb.' * rb
        gA(:, b) = Jb.' * rb;                     % [nMat x 1]

        % Gauss-Newton Hessian term H1
        wgn = (lambda_true2(:,b) + lambda_contrast(:,b)) ./ (lambda_contrast(:,b).^2); % [nBins x 1]
        H1 = Jb.' * bsxfun(@times, Jb, wgn);      % [nMat x nMat]

        % curvature term H2
        c = (rb.' * w);                           % [1 x nE]
        t = bsxfun(@times, c, tb);                % [1 x nE]  (c .* expb .* dE)
        Ft = bsxfun(@times, F, t);                % [nMat x nE]
        H2 = Ft * F.';                            % [nMat x nMat]

        % total Hessian for batch b (symmetry guard)
        Hb = H1 + H2;
        HA(:,:,b) = 0.5 * (Hb + Hb.');
    end
end