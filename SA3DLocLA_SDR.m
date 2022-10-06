function [uRfn,uSDR] = SA3DLocLA_SDR(theta,gamma,S,Q)
% [uRfn,uSDR] = SA3DLocLA_SDR(theta,gamma,S,Q)
% 
% Solution for 3-D localization in arbitrary geometry using space angles 
% from linear arrays, by SDR with refinement.
%
% Input parameters:
% theta:    (M x 1), space angle (SA) measurements, in radian.
% gamma:    (3 x M), directions of linear arrays.
% S:        (3 x M), positions of linear arrays.
% Q:        (M x M), measurement noise covariance matrix.
%
% Output parameters:
% uSDR:      (3 x 1), source position estimate by SDR.
% uRfn:      (3 x 1), refined source positon estimate.
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

[N,M] = size(S);

h = diag(gamma'*S);
G = [gamma', -diag(cos(theta))];
B = eye(M);

for itr = 1:2
    W = inv(B*Q*B);
    A = [G'*W*G, -G'*W*h;
        -h'*W*G, h'*W*h];

    cvx_begin quiet
    cvx_solver sdpt3    % TSP 2022
        variable P(M+N,M+N) symmetric
        variable p(M+N,1)
        SM = [P,p;p',1];
        minimize trace(SM*A)
        subject to
            for i = 1:M
                P(i+N,i+N) == trace(P(1:N,1:N)) - 2*S(:,i)'*p(1:N) + S(:,i)'*S(:,i);
                p(i+N) >= 0;
            end
            for i=1:M-1
                for j=i+1:M
                    P(i+N,j+N) >= abs(trace(P(1:N,1:N)) - S(:,i)'*p(1:N) - S(:,j)'*p(1:N) + S(:,i)'*S(:,j));
                end
            end
            SM == semidefinite(M+N+1);
    cvx_end

    u2 = sqrt(diag(SM(1:N,1:N)));
    uSDR = func_findu(theta,gamma,S,Q,u2,M);

    r = sqrt(sum((S-uSDR).^2,1))';
    B = -diag(r.*sin(theta));
end


% refinement stage
uRfn = uSDR;

for iter = 1:2
    r = sqrt(sum((S-uRfn).^2,1))';
    g = h - gamma'*uRfn + diag(cos(theta))*r;

    Ht(1:N,1:N) = eye(N);
    for i = 1:M
        Ht(N+i,:) = (uRfn-S(:,i))'/norm(uRfn-S(:,i),2);
    end
    H2 = G*Ht;
    du2 = -(H2'*W*H2)\H2'*W*g;
    uRfn = uRfn - du2;
    
    r = sqrt(sum((S-uSDR).^2,1))';
    B = -diag(r.*sin(theta));
    W = inv(B*Q*B);
end

%varargout = {pos,ue};
end




function uu = func_findu(theta,gamma,S,Q,pos,M)
    sgn = [1     1     1
           1     1    -1
           1    -1     1
           1    -1    -1
          -1     1     1
          -1     1    -1
          -1    -1     1
          -1    -1    -1]';
    for j = 1:8
        u = sgn(:,j).*pos;
        for i = 1:M
            est_theta(i,1) = acos(gamma(:,i)'*(u-S(:,i))/norm(u-S(:,i)));
        end
        J(j) = (theta - est_theta)'/Q*(theta - est_theta);
    end

    [~,ind] = min(J);
    uu = sgn(:,ind).*pos;
end