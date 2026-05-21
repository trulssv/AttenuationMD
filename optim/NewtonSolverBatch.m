function [A_opt, out] = NewtonSolverBatch(AB, fB, AR, fR, F, AC0,  weights,dE, opts)
% minimize_KL_Newton: Minimize D_JS(Pm, Q(A)) where Q(A) has Poisson rates
% tilde_lambda_i(A) = ∫ w_i(E) * exp(-sum_j A_j f_j(E)) dE  (discretized).
%
%   Inputs:
%   AB          : [nBatch x nMat] Residual
%   AR          : [nBatch x nMat] Background
%   A_init      : [nBatch x nMat] initial guess for A (nonnegative)
%   weights     : [nBins x nE] w_i(E)
%   F           : [N x nE] basis f_j(E)
%   dE          : scalar or [1 x nE] energy widths
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
%   A_opt : [nBatch x 1] optimized A
%   out   : struct with diagnostics (history of f, grad norms, etc.)

    % --- robust opts initialization ---
    if nargin < 10 || isempty(opts)
        opts = struct();
    end

    function val = local_getopt(name, def)
        if isfield(opts, name) && ~isempty(opts.(name))
            val = opts.(name);
        else
            val = def;
        end
    end

    maxIter   = local_getopt('maxIter',   300);
    gtol      = local_getopt('gtol',      1e-2);
    xtol      = local_getopt('xtol',      1e-8);
    ftol      = local_getopt('ftol',      1e-4);
    mu0        = local_getopt('mu0',       1e-3);
    floorRate = local_getopt('floorRate', 1e-12);
    verbose   = local_getopt('verbose',   0);
    mu_reg    = local_getopt('mu_reg',     1e-8);


    % Ensure shapes
    [nBins, nE] = size(weights);
    [nMat, nE2]= size(F);
    [nBatch, nMat2] = size(AC0);
    if nE ~= nE2, error('weights and F must have same # of energy bins'); end
    if nMat ~= nMat2, error('Initial guess and F must have the same # of materials'); end

    if or(size(AB)~=size(AR), size(AB)~=size(AC0)), error('AB, AR, AC0 must have the same size'); end
    
    %% Reshape data
    
    %  max x mat x bin x E  x batch
    
    weights = reshape(weights, [1, 1, nBins, nE]);
    A = reshape(AC0', [1, nMat, 1, 1, nBatch]);

    
    %% Compute attenuation, updated wieghts and lambda_true
    
    attenuationB0 = exp(-AB* fB)'; % Attenuation from AB: nE x nBatch 
    
    attenuation_B(1, 1, 1, :, :) = attenuationB0; %1 x 1 x 1 x nE x nBatch
    
    
    attenuationR0 = exp(-AR* fR)'; % Attenuation from AR: nE x nBatch
    attenuation_R(1, 1, 1, :, :) = attenuationR0; %1 x 1 x 1x nE x nBatch

    
    weights_updated = attenuation_B .* weights;             % weights combined with attenation from AB: 1x 1 x nBins x nE x nBatch
    integrand_lambda_true = weights_updated.*attenuation_R * dE;  % 1 x 1 x nBins x nE x nBatch
    lambda_true = sum(integrand_lambda_true, 4);            % lambdas according to CatSim: 1x 1 x nBins x 1 x nBatch
    
    
    
    %% Define variables
    
    step_norm = zeros(1, 1, 1, 1, nBatch);
    
    %% Initial regularization 
    
    A_batch = max(A, floorRate);
    [~, gA_batch, HA_batch] = NewtonStepBatch(A_batch, lambda_true, weights_updated, F, dE, nMat, nBins, nBatch, floorRate);
            
    Hscale_batch = max(1, max(abs(HA_batch.*repmat(eye(nMat), 1, 1, 1, 1, nBatch)), [], [1, 2]));
    mu_batch = mu0 * Hscale_batch;
    
    for i=1:maxIter
        
        
        
        [f_batch, gA_batch, HA_batch] = NewtonStepBatch(A_batch, lambda_true, weights_updated, F, dE, nMat, nBins, nBatch, floorRate);
        
        f_trial_batch = zeros(size(f_batch));
        s_batch = zeros(1,  nMat, 1, 1, nBatch);
        
        Hscale_batch = max(1, max(abs(HA_batch.*repmat(eye(nMat), 1, 1, 1, 1, nBatch)), [], [1, 2]));
        reg_batch = max(mu_batch + mu_reg, 1e-12 * Hscale_batch);
        Hreg_batch = HA_batch + reg_batch .* eye(nMat);
        
        %% regularize Hessian
       
        
        
        %% stop if all small gradients
        
        if mean(vecnorm(gA_batch), 'all') < gtol / 10
            if verbose, fprintf('BatchNewton converged in %d iters\n', i-1); end
        break;
        end
                
        for k=1:nBatch
            A = reshape(A_batch(:, :, :, :, k), [nMat, 1]);
            gA = reshape(gA_batch(:,:, :, :, k), [nMat, 1]);
            HA = HA_batch(:, :, :, :, k);
            f = f_batch(:, :, :, :, k);
            Hreg = Hreg_batch(:, :, :, :, k);
            reg =reg_batch(:, :, :, :, k);
            
             % Compute regularization for THIS batch element with current mu
            
            
            
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
                    Hreg = HA + jitter * eye(nMat);
                    attempts = attempts + 1;
                end
            end
            if attempts == 10
                warning('Could not make Hessian SPD after jitter attempts; using backslash as fallback.');
                p = - Hreg \ gA;  % fallback
            end
            
            %% Armijo backtracking in B-space
        
            c1 = 1e-4; % Armijo paramter
            [alpha, f_trial] = ArmijoBacktracking(); % Armijo etc using f(B)->f(A=exp(B))
            s = alpha * p;            
            
            %% Accept step
            
        
            A_batch(:, :, :, :, k) = A + s;
            f_batch(k) = f;
            f_trial_batch(:, :, :, :, k) = f_trial;
            s_batch(:, :, :, :, k) = s;
            f = f_trial;           % update current objective to accepted trial

            step_norm(:, :, :, :, k) = norm(s);
            
            
        end
        
        %% Update damping with vectorized Levenberg-Marquardt rule
        
        LevenbergMarquardtUpdate();
        
        %% stopping on step size
        if mean(step_norm, 'all') < xtol
            if verbose, fprintf('Converged: step small\n'); end
            break;
        end
        
         %% stopping on step size
        if mean(f_batch) < ftol
            if verbose, fprintf('Converged: f small\n'); end
            break;
        end
        
        if verbose
            fprintf('iter: %d mean(gA)= % 0.6e step_norm: % 0.6e, mean(f)= % 0.6e\n', i-1, mean(vecnorm(gA_batch)), mean(step_norm), mean(f_batch))
        end
        
    end
    
    out = squeeze(f_batch);
    A_opt = squeeze(A_batch);

function [alpha, f_trial] = ArmijoBacktracking()
    
    alpha = 1;
    f_curr = f;
    A_curr = A;
    gA_curr = gA;
    lambda_curr = squeeze(lambda_true(:, :, :, :, k));
    weights_curr = squeeze(weights_updated(:, :, :, :, k));
    maxbacktracks = 20;
   
    for l=1:maxbacktracks
        A_trial = A_curr + alpha * p;
        
        [f_trial,~, ~]  = NewtonStep(A_trial, lambda_curr, weights_curr , F, dE, floorRate);
        
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
    
    %% Reshaping
    
    s_batch = reshape(s_batch, [nMat, 1, 1, 1, nBatch]);
    
    %% Compute step quality
    
    s = p * alpha; % Acutal step
    
    actual_reduction = f_batch - f_trial_batch;
    predicted_reduction = - pagemtimes(gA_batch, s_batch) - 0.5 * pagemtimes(pagemtimes(s_batch, 'transpose', HA_batch, 'none'), s_batch);
    
    eps = 1e-10;
    predicted_reduction = max(eps, predicted_reduction);
    
    %% Calculate rho
    
    rho_batch = actual_reduction ./ predicted_reduction;
    
    %% Update mu based on rho
    
    shrinkMask = rho_batch > rho_up;
    growMask   = rho_batch < rho_down;

    mu_batch = max(mu_min, mu_batch .* (factor_up.^shrinkMask));
    mu_batch = min(mu_max, mu_batch .* (factor_down.^growMask));
    
end

end

function [f, gA, HA] = NewtonStepBatch(A, lambda_true, weights, F, dE, nMat, nBins, nBatch, floorRate)
    % Compute value, gradient (w.r.t A), and Hessian (w.r.t A) of
    
        
    Attenuation_C = exp(-pagemtimes(F, 'transpose', A, 'transpose'));     %Contrast attenuation: [nE x 1 x 1 x 1 x nBatch]
    
    Attenuation_C = permute(Attenuation_C, [2, 3, 4, 1, 5]); % 1 x 1 x nE x nBatch
  
   
   % A should be [nMat x 1 x nBatch] for multiplication with F [nMat x nE]
%    A_reshaped = permute(A, [2, 4, 5, 1, 3]); % [nMat x 1 x nBatch]
%    F_expanded = repmat(F, [1, 1, nBatch]); % [nMat x nE x nBatch]
%    z = pagemtimes(A_reshaped, 'transpose', F_expanded, 'none'); % [1 x nE x nBatch]
%    Attenuation_C = exp(-z);
   
%    
%    fprintf('F shape: [%s]\n', num2str(size(F)));
%    fprintf('A shape: [%s]\n', num2str(size(A)));
%    fprintf('Attenuation_C shape: [%s]\n', num2str(size(Attenuation_C)));
    
    F = reshape(F, [1, size(F, 1), 1, size(F, 2)]);      %1 x nMat x 1 x nE x 1
    
    
    
    % Discretized integrals
    S = weights .* Attenuation_C * dE;   % [1 x nBins x nE x nBatch]

    % Rates
    lambda_contrast = sum(S, 4);   % [1 x nBins x 1 x nBatch]
    lambda_contrast = max(lambda_contrast, floorRate);      % positivity floor
    
    %% Minimize symmetrised KL-divergece: KL(P || Q) + KL(Q || P).
    
    
    % Jacobian of lambda_contrast w.r.t A:
    
    Jint = - F.*S;                              % nMat x nMat x nBins x nE xnBatch 
    J = sum(Jint, 4);                           % Integrate over energy: nMat x nMat x nBins x 1 x nBatch
    
    % Hessian terms
     
    F2  = permute(F, [2, 1, 3, 4]);
    DJint = sum(-Jint .*F2, 4);
        
    % Residual weight vector
    
    r = 1 - (lambda_true ./ lambda_contrast) + log(lambda_contrast ./ lambda_true);   % [1 x 1 x nBins x 1 x nBatch]

    % Calculate batch divergence
    
    f = sum( lambda_true .* (log(lambda_true ./ lambda_contrast)) + lambda_contrast .* (log(lambda_contrast ./ lambda_true)), 3); % 1 x 1 x 1 x 1 x nBatch 
    f = abs(f); % Ensure f is positive for numerical stability.

    % Hessian w.r.t A:

    h1 = (lambda_true + lambda_contrast) ./ (lambda_contrast.^2);           % [1 x 1 x nBins x 1 x nBatch]

    % Gradient w.r.t A: gA = J' * r
    gA = sum(J.* r, 3);                                     % [nMat x nBin x 1 x nBatch]
    
    
    % 1st Hessian term 
    
    Jh1 = J.*h1;
    
    H1 = pagemtimes(permute(Jh1, [2, 3, 1, 4, 5]), permute(J, [3, 2, 1, 4, 5]));
        
    % 2nd Hessian term
    
    H2 = sum(DJint.*r, 3);

    HA = (H1 + H2);
    
    % Numerical guarding (optional, keeps symmetry)
    
    HA = 0.5*(HA + permute(HA, [2, 1, 3, 4, 5]));
    
     
end

function [f, gA, HA] = NewtonStep(A, lambda_true, weights, F, dE, floorRate)
    % Compute value, gradient (w.r.t A), and Hessian (w.r.t A) of
    
    % Shapes
    [M, ~] = size(weights);
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
    
    % Residual weight vector
    r = 1 - (lambda_true ./ lambda_model) + log(lambda_model ./ lambda_true);           % [M x 1]

    % Value
    f = sum( lambda_true .* (log(lambda_true ./ lambda_model)) + lambda_model .* (log(lambda_model ./ lambda_true)));
    f = abs(f); % Ensure f is positive for numerical stability.


    % Hessian w.r.t A:
    % GN term: J' * diag(w) * J with w_i = lambda_i / tilde_lambda_i^2
    wgn = (lambda_true + lambda_model) ./ (lambda_model.^2);           % [M x 1]

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