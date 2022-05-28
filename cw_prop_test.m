clear; clc;

% parameters
mu = 398600;            % [km3/s2]
rRef = 6600;            % 222km altitude [km]
n = sqrt(mu/rRef^3);

deltaT = 0.01;
% t = 0.5*3600;
t = 2.75;

r0 = [0 0 0 -0.0656 -0.0002 0]';

x_hist = cw_propagator(r0,n,t,deltaT);

