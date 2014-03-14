% solveDDGPS.m	(actual file name: solveDDGPS.m)
%
% this CDGPS utility function computes delta X from the equation: 
% l_k = A_k*x + lambdaL1*bias
%
% input:    'biases'            the biases calculated. this could be
%                               real-value biases or the integer ambiguity 
%                               lattice points. [N-1] by 1
%
%           'A'                 A matrix calculated in 'calculateQtilda.m'
%                               that is used in the above equation.
%
%           'l_k'               l_k matrix calculated in 'l_kmatrix.m' used
%                               for the solving of del_x
%
% output:   'del_x'             3-1 matrix containing delta X, Y, and Z of
%                               the ECEF calculations
%
function del_x = solveDDGPS(biases,A, l_k)

% download the constant file
    constant;
% create Q_dd
    Q_dd = 2.*ones(length(A)) + 2.*eye(length(A));
% calculate del_x
    del_x = inv((A'*inv(Q_dd)*A))*A'*inv(Q_dd)*(l_k - lambdaL1*biases);  
end