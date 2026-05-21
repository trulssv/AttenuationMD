function [A_opt, out] = MDSolve(fB, AB, fR, AR, fC, AC0, ...
    weights, eActVectorkeV, mode)
    
    
    assert(nargin >=6, 'Error: not enough arguments')
    
    if nargin < 8
        % Default
        
        Emin = 15;
        Emax = 120;
        dE = 1;
        eActVectorkeV  = struct('E', {Emin:dE:Emax}, 'dE',{dE});
        
        weights = GetWeightfunctions(eActVectorkeV.E);
        
        mode = 'simple';
        
        
    elseif nargin == 8
        mode = 'simple';
    end
    

    
    F             = fC;             % [N x nE]
    dE            = eActVectorkeV.dE;     % scalar or 1 x nE
    A_init        = AC0;        % [N x 1] initial guess for A^C
    
    %% Solver hyperparams
    
    opts = struct('maxIter',100,'gtol',1,'xtol',1e-6,'ftol',1e-4, 'floorRate', 1e-8, 'mu0',1e-3, 'mu_reg', 1e-7, 'verbose',0);    
    
    switch mode
        
        case 'reparameterized'
            
            weights_model = GetUpdatedWeights(weights, AB, fB); 
            lambda_true     = GetLambda(weights_model, eActVectorkeV, fR, AR);
            [A_opt, out] = NewtonSolverReparameterization(lambda_true, weights_model, F, dE, A_init, opts);

        case 'vectorized'
            
            %% Preprocess, identify non-zero values. 
            
            ep = 1e-2;
            
            nonzeros_idx = sum(abs(AR), 2) > ep;
            ARnz = AR(nonzeros_idx, :);
            ABnz = AB(nonzeros_idx, :);
            AC0nz = AC0(nonzeros_idx, :);
            
            [Ptot, ~] = size(AC0);
            
            [P, nMat] = size(AC0nz);
            
            batch_size = 1000;
            num_batches = ceil(P / batch_size);
                        
            % Execute batches in parallell
            
            A_result_batches = cell(1, num_batches);
            
            parfor b=1:num_batches
               
                fprintf('Processing batch %d of %d total batches \n', b, num_batches)
                
               idx = (b-1)*batch_size + (1:min(batch_size, P-(b-1)*batch_size));
               
               
               AB_batch = ABnz(idx, :);
               AR_batch = ARnz(idx, :);
               AC0_batch = AC0nz(idx, :);
               
               
               [A_opt_batch, ~] = NewtonSolverBatch(AB_batch, fB, AR_batch, fR, F, AC0_batch,  weights,dE, opts); % vectorized energy sums + small local loops
               A_result_batches{b} = {A_opt_batch, idx};
       
               
            end
            
            % stitch back together
            
            A_result = zeros(P, nMat);
            for b = 1:num_batches
                idx = A_result_batches{b}{2};
                A_result(idx, :) = squeeze(A_result_batches{b}{1})';
            end
            
            out = NaN;
            
            A_opt = zeros(Ptot, nMat);
            A_opt(nonzeros_idx, :) = A_result;
            
            
            
        case 'simple'
            weights_model = GetUpdatedWeights(weights, AB, fB); 
            lambda_true     = GetLambda(weights_model, eActVectorkeV, fR, AR);
            [A_opt, out] = NewtonSolver(lambda_true, weights_model, F, dE, A_init, opts);
        case 'single'
                    weights_model = GetUpdatedWeights(weights, AB, fB); 
                    lambda_true     = GetLambda(weights_model, eActVectorkeV, fR, AR);
                    [A_opt, out] = NewtonSolver(lambda_true, weights_model, F, dE, A_init, opts);
                    out = out.final_f;
    end
    
end



function updated_weights = GetUpdatedWeights(weigths, AB, fB)
    attenuation = exp(-AB * fB);
    updated_weights = weigths .* attenuation;
end

function lambda = GetLambda(weigths, eActVectorkeV, fR, AR)
     attenuation = exp(-AR * fR);
     integrand = weigths .* attenuation * eActVectorkeV.dE;
     lambda = sum(integrand, 2);
end
