function [CRB] = SA3DLocLACRLB(u,gamma,S,Q)
% [CRB] = SA3DLocLACRLB(u,gamma,S,Q)
%
% CRLB for source position using space angle measurements.
%
% Input parameters:
% u:        (3 x 1), source position.   
% gamma:    (3 x M), directions of linear arrays.
% S:        (3 x M), positions of linear arrays.
% Q:        (M x M), measurement noise covariance matrix.
%
% Output parameters:
% CRB:      (3 x 3), CRLB matrix of source position estimate
%
% Reference:
% Y. Sun, K. C. Ho, L. Gao, J. Zou, Y. Yang, and L. Chen, "Three 
% dimensional source localization using arrival angles from linear arrays: 
% analytical investigation and optimal solution," IEEE Trans. Signal 
% Process., vol. 70, pp. 1864-1879, 2022.
%                                                                        
% Yimao Sun and K. C. Ho   05-08-2022
%
%       Copyright (C) 2022
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA.
%       hod@missouri.edu
%

[K,M] = size(S);             % K=dimension, S=number of sensors
r = sqrt(sum((u-S).^2,1))';

for i = 1:M
    theta(i,1) = acos(gamma(:,i)'*(u-S(:,i))/norm(u-S(:,i)));
            A(i,:) = gamma(:,i)'*((u-S(:,i))*(u-S(:,i))'-r(i)^2*eye(K)) / (r(i)^3*abs(sin(theta(i))));
end

CRB = inv(A'/Q*A);
