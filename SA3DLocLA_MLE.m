function uMLE = SA3DLocLA_MLE(theta,gamma,S,Q,uInit,varargin)
% uMLE = SA3DLocLA_MLE(theta,gamma,uInit,S,Q,varargin)
%
% Maximum Likelihood Estimator (MLE) for source localization in 3-D 
% using space angles from linear arrays.
%
% Input parameters:
% theta:    (M x 1), space angle (SA) measurement vector, in radian.
% gamma:    (3 x M), directions of linear arrays.
% S:        (3 x M), positions of linear arrays.
% Q:        (M x M), measurement noise covariance matrix.
% uInit:    (3 x 1), initial solution guess.
% varargin: number of iterations if not empty.
%
% Output parameter:
% uMLE:     (3 x 1), source localization estimate
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
            
warning off;

[~,M] = size(S);  % M is the number of sensors

if isempty(varargin)
    N = 20;      % default max.             number of iterations
else
    N = varargin{1};
end

for k = 1:N
    
    for i = 1:M
        dp = uInit - S(:,i);
        r = norm(dp);
        l = gamma(:,i)'*dp;
        h = sqrt(r^2 - l^2);
        K(i,:) = gamma(:,i)'/h - l*dp'/r^2/h;
        thetaHat(i,1) = acos(gamma(:,i)'*dp/r);
    end
    
    delta = -(K'/Q*K)\K'/Q*(theta-thetaHat);

    if norm(delta) < 0.0001 || sum(isnan(delta))
        break;
    else
        uInit = uInit + delta;
    end
end

uMLE = uInit;

end