function [x, iteration_counter] = NR_function(F, J, x, eps)
% Solves nonlinear system F = 0 by Newton-Raphson's method

% INPUT:    F   ... nonlinear system of equations as F = 0
%           J   ... Jacobian matrix of the system F
%           x   ... starting value of x
%           eps ... tolerance of the norm of F
% OUTPUT:   x   ... x as result fullfiling condition norm of F < eps
%           iteration_counter ...  number of iterations of NR method to get
%           the result

F_value = F(x);
F_norm = norm(F_value);
iteration_counter = 0;

while F_norm > eps && iteration_counter < 100
    delta = J(x)\-F_value;
    x = x + delta;
    F_value = F(x);
    F_norm = norm(F_value);
    iteration_counter = iteration_counter + 1;
end
if F_norm > eps
    iteration_counter = -1;
end

end

