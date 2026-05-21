function M = maskedLowPass(Mk, mask, kernelsize)
    h = ones(kernelsize, kernelsize);
    num = conv2(Mk, h, 'same');
    denom = conv2(mask, h, 'same');
    denom(denom==0) = NaN;
    M = num./denom;
    M(~mask) = 0;

end