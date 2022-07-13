% function sse = fit_singleGaussian(params, Input, Actual_Output)
%   Fitted_Curve = A * exp( - (Input - Mu).^2 / SD.^2 );

function sse = fit_singleGaussian(params, Input, Actual_Output)
    A = params(1);
    Mu = params(2);
    SD = params(3);
    
    Fitted_Curve = A * exp( - (Input - Mu).^2 / SD.^2 );
    
    Error_Vector = Fitted_Curve - Actual_Output;
    
    % typically minimize the sum of square error of curve fitting as:
    sse=sum(Error_Vector.^2);
