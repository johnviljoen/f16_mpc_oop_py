close all
clear all
clc

npos = 0;
epos = 0;
h = 10000;
phi = 0;
theta = 0;
psi = 0;
V = 700;
alpha = 0;
beta = 0;
p = 0;
q = 0;
r = 0;

% npos_dot = 0;
% epos_dot = 0;
% h_dot = 0;
% phi_dot = 0;
% theta_dot = 0;
% psi_dot = 0;
% V_dot = 0;
% alpha_dot = 0;
% beta_dot = 0;
% p_dot = 0;
% q_dot = 0;
% r_dot = 0;

T = 2886.6468;
dh = -2.0385;
da = -0.087577;
dr = -0.03877;
lef = 0.3986;
fi_flag = 1;

x = [npos;epos;h;phi;theta;psi;V;alpha;beta;p;q;r;T;dh;da;dr;lef;fi_flag];

xdot = zeros(18,1);

xdot = nlplant(x)