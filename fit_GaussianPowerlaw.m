function sse = fit_GaussianPowerlaw(params, Input, Actual_Output)
    A = params(1);
    Mu = params(2);
    SD = params(3);
    B = params(4);
    k = params(5);
    x0 = params(6);
    
    Fitted_Curve = A * exp( - (Input - Mu).^2 / SD.^2 ) + B * Input.^(-k) .* (Input > x0);
    
    Error_Vector = Fitted_Curve - Actual_Output;
    
    % typically minimize the sum of square error of curve fitting as:
    sse=sum(Error_Vector.^2);
