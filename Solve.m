function AC_est = Solve(AB, AR)
 

    sz = size(AB{1});

    Emin = 15;
    Emax = 120;
    dE = 1;
    eActVectorkeV  = struct('E', {Emin:dE:Emax}, 'dE',{dE});
    
    load("weights.mat", "weights")
    load('MuGdimm.mat', 'MuGdimm')
    load('MuIimm.mat', 'MuIimm')
    load('MuWaimm.mat', 'MuWaimm')

    % weights = GetWeightfunctions(eActVectorkeV.E);
    % 
    % MuPEimm = 0.1*GetMu('polyethylene', eActVectorkeV.E);
    % MuPVCimm = 0.1*GetMu('pvc', eActVectorkeV.E);
    % MuIimm = 0.1*GetMu('I', eActVectorkeV.E);
    % MuGdimm = 0.1*GetMu('Gd', eActVectorkeV.E);

    fB = [MuWaimm; MuIimm];
    fR = [MuWaimm; MuIimm];
    fC = [MuIimm; MuGdimm];
    
    AB = [ AB{1}(:), AB{2}(:)];
    AR = [AR{1}(:), AR{2}(:)];
    AC0 = 0.05 * ones(size(AB));
    
    tic
    [AC, ~] = MDSolve(fB, AB, fR, AR, fC, AC0, weights, eActVectorkeV, 'vectorized');
    toc
    
    Ipred_sino = AC(:, 1);
    Gdpred_sino = AC(:, 2);
    
    % Post-processing:
    
    % remove nans
    
    %Ipred_sino(isnan(Ipred_sino)) = 0;
    %Gdpred_sino(isnan(Gdpred_sino)) = 0;
    
    
    % remove negative values
    
    %minthresholdI = -0.1;
    %minthresholdGd = -0.01;
    
    %Ipred_sino(Ipred_sino < minthresholdI) = 0;
    %Gdpred_sino(Gdpred_sino < minthresholdGd) = 0;
    
    % max threshold
    
    %maxthresholdI = 0.15;
    %maxthresholdGd = 0.03;

    %Ipred_sino(Ipred_sino > maxthresholdI) = 0;
    %Gdpred_sino(Gdpred_sino > maxthresholdGd) = 0;
    
    % Reconstruct as sinos
    
    Ipred_sino = reshape(Ipred_sino, sz);
    Gdpred_sino = reshape(Gdpred_sino, sz);
    
    AC_est = {Ipred_sino, Gdpred_sino};


    
end