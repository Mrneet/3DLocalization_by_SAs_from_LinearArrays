% This code generates Figs. 10-21 in the paper
% Y. Sun, K. C. Ho, L. Gao, J. Zou, Y. Yang, and L. Chen, "Three
% dimensional source localization using arrival angles from linear arrays:
% analytical investigation and optimal solution," IEEE Trans. Signal
% Process., vol. 70, pp. 1864-1879, 2022.
%
% Note: the code has 20 ensemble runs only to reduce computation time.
%       Please change the variable "mon" to 500 for obtaining the results
%       in the paper.
%   
% Yimao Sun and K. C. Ho   05-08-2022
%
%       Copyright (C) 2022
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA.
%       hod@missouri.edu
%

clear; %close; clc;
warning off;

rng('default');

clor = [0, 114, 189; 217, 83, 25; 237, 177, 32; 126, 47, 142; 119, 172, 48; 77, 190, 238; 162, 20, 47]/256;

u = [424;519;375];  % true source position

% -- scenarios to generate varous figures in the paper --
scen = 'general';   % Figs. 10-12
% scen = 'coplanar';  % Figs. 13-15
% scen = 'uniDirec';  % Figs. 16-18
% scen = 'uniCopla';  % Figs. 19-21
% scen = 'uniDirec2';
% scen = 'uniCopla2';

% -- settings --
mon=20;     % 500               % number of ensemble runs


% ------------------------------------------------------
switch scen
    case 'general'
        M = 7;          % number of sensors
        aTmp = ((0:1/M:0.99))'*2*pi;
        S = [300, 187,  -67, -270, -270,  -67,  187;
            0,   235,  292,  130, -130, -292, -235;
            94,  122, -112,  124,   39, -121,  -66]; % sensor positions
        alpha = aTmp;   % azimuth of array direction
        beta = [0;0.813;0.228;2.460;2.270;0.438;1.500]; % elevation of array direction
        nsepwr_dB = -70:10:10; % noise power in log scale (10 log10(nsePwr))
        plotLabel = 1:3;
        ylimMSE = [0.05,1e8];
        ylimBias = [0.1,1e4];
    case 'coplanar'
        M = 7;
        aTmp = ((0:1/M:0.99))'*2*pi;
        S = [ 236,    3,   37,  173,  -89, -167,  193;
            -254, -250, -114,  213,  -15,  105,  -23;
            0,       0,    0,    0,    0,    0,    0];
        alpha = aTmp;
        beta = zeros(M,1);
        nsepwr_dB = -70:10:10;
        plotLabel = 1:3;    % 1:6;
        ylimMSE = [0.05,1e8];
        ylimBias = [0.1,1e4];
    case 'uniDirec'
        M = 8;
        aTmp = 2*pi/M*(1:M)'-pi/M;
        S = [277, 115, -115, -277, -277, -115,  115,  277;
            115, 277,  277,  115, -115, -277, -277, -115;
            94,  137,  -23,   54,  -67,  -18,   63,  140];
        alpha = zeros(M,1);
        beta = zeros(M,1);
        nsepwr_dB = -80:10:0;
        plotLabel = 1:3;    % [1:5,7];
        ylimMSE = [0.05,1e8];
        ylimBias = [0.1,1e4];
    case 'uniCopla'
        M = 8;
        aTmp = 2*pi/M*(1:M)'-pi/M;
        S = [277, 115, -115, -277, -277, -115,  115,  277;
            115, 277,  277,  115, -115, -277, -277, -115;
            0,     0,    0,    0,    0,    0,    0,    0];
        alpha = zeros(M,1);
        beta = zeros(M,1);
        nsepwr_dB = -80:10:0;
        plotLabel = 1:3;    % [1:5,7];
        ylimMSE = [0.01,1e8];
        ylimBias = [0.1,1e4];
    case 'uniDirec2'
        M = 8;
        aTmp = 2*pi/M*(1:M)'-pi/M;
        S = [277, 115, -115, -277, -277, -115,  115,  277;
            115, 277,  277,  115, -115, -277, -277, -115;
            94,  137,  -23,   54,  -67,  -18,   63,  140];
        alpha = rand*pi*ones(M,1);
        beta = zeros(M,1);
        nsepwr_dB = -80:10:0;
        plotLabel = 1:3;    % [1:5,7];
        ylimMSE = [0.01,1e8];
        ylimBias = [0.1,1e4];
    case 'uniCopla2'
        M = 8;
        aTmp = 2*pi/M*(1:M)'-pi/M;
        S = [277, 115, -115, -277, -277, -115,  115,  277;
            115, 277,  277,  115, -115, -277, -277, -115;
            0,     0,    0,    0,    0,    0,    0,    0];
        alpha = rand*pi*ones(M,1);
        beta = zeros(M,1);
        nsepwr_dB = -80:10:0;
        plotLabel = 1:3;    % [1:5,7];
        ylimMSE = [0.01,1e8];
        ylimBias = [0.1,1e4];
end

gamma = [cos(alpha).*cos(beta),sin(alpha).*cos(beta),sin(beta)]';   % array direction

% *****************************************
% localization geometry
figure;
plot3(S(1,:),S(2,:),S(3,:),'o','linewidth',1.5); hold on;grid on;
plot3(u(1),u(2),u(3),'x','linewidth',1.5); axis equal;
rr = 30;
s_start = S - rr*gamma;
s_end = S + rr*gamma;
for m = 1:M
    plot3([s_start(1,m),s_end(1,m)],[s_start(2,m),s_end(2,m)],[s_start(3,m),s_end(3,m)],'k-','linewidth',1.5);
end
legend('Arrays','Source','Attitude','Location','best','FontSize',11);
xlabel('x(m)','fontsize',12);ylabel('y(m)','fontsize',12);zlabel('z(m)','fontsize',12);
% *****************************************

% true space angles
for i = 1:M
    thetao(i,1) = acos(gamma(:,i)'*(u-S(:,i))/norm(u-S(:,i)));
end

nsepwr = 10.^(nsepwr_dB/10);    % noise power

for m = 1:mon                   % pre-generate noise pattern
    nseTmp(:,m) = randn(M,1);
end
nse = nseTmp - mean(nseTmp,2);

atmp = rand(M,1)+1/9;
Ra = diag(roundn(atmp/mean(atmp),-2));

for n = 1:length(nsepwr)        % for each noise level
    disp([num2str(n),'/',num2str(length(nsepwr))]);
    Q = nsepwr(n)*Ra;
    CRB = SA3DLocLACRLB(u,gamma,S,Q);
    crlb(n) = trace(CRB);

    rng('default');
    for m = 1:mon               % ensemble runs
        theta = thetao + sqrtm(Q)*nse(:,m);
        ia = 1;

        % Proposed (SDR + refinement) solution
        [ uRE, uSDR] = SA3DLocLA_SDR( theta,gamma,S,Q );
        uEst(:,ia,m) = uSDR; ia = ia + 1;
        uEst(:,ia,m) = uRE; ia = ia + 1;

        % MLE method initialized by SDR solution
        uGN = SA3DLocLA_MLE(theta,gamma,S,Q,u);
        %uGN = AOA3DLocLA_MLE(theta,gamma,S,Q,uSDR);
        uEst(:,ia,m) = uGN; ia = ia + 1;
    end
    mse(:,n) = mean(sum((uEst-u).^2,1),3);
    bias(:,n) = mean(sqrt(sum((uEst-u).^2,1)),3);
end

symbs = ['o','^','x','s','*','v','+'];
names = {'SDR-Ini','Refined','MLE','SDR-Zou','ADMM','WLS-Zou','WLS-Wei'};

% -- plot results --
figure;
for i = plotLabel
    semilogy(nsepwr_dB,(mse(i,:)),symbs(i),'Linewidth',1.5,'color',clor(i,:),'DisplayName',names{i}); hold on; grid on;
end
semilogy(nsepwr_dB,(crlb),'-','Linewidth',1.5,'DisplayName','CRLB');
legend('Show','Location','Northwest','fontsize',11);
xlabel('10log(\sigma^2(rad^2))','fontsize',12); ylabel('MSE(u) (m^2)','fontsize',12);
ylim(ylimMSE);

figure;
for i = plotLabel
    semilogy(nsepwr_dB,(bias(i,:)),symbs(i),'Linewidth',1.5,'color',clor(i,:),'DisplayName',names{i}); hold on; grid on;
end
legend('Show','Location','Northwest','fontsize',11);
xlabel('10log(\sigma^2(rad^2))','fontsize',12); ylabel('Bias(u) (m)','fontsize',12);
ylim(ylimBias);
