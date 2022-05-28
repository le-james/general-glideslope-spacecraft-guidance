clear; clc;

% choosen by mission requirements
r0 = [500; 0; 200];       % chaser inital position relative to the target
rT = [100; 0; 0];    % goal position of the chaser relative to the target
rho_0 = norm(r0-rT);    % magnitude of r0     
rho_dot_0 = -1;       % inital commanded velocity [m/s]
rho_dot_T = -0.3;       % velocity to arrive at target position [m/s]

% compute unit vector components
x0 = r0(1);
y0 = r0(2);
z0 = r0(3);
xT = rT(1);
yT = rT(2);
zT = rT(3);
cos_alpha = (x0-xT)/rho_0;
cos_beta = (y0-yT)/rho_0;
cos_gamma = (z0-zT)/rho_0;
% compute unit vector components
rho_unit_vec = [cos_alpha cos_beta cos_gamma]'; % rho unit direction vector


a = (rho_dot_0-rho_dot_T)/rho_0;        % a<0

T = (1/a)*log(rho_dot_T/rho_dot_0);     % total time of rendezvous trajectory
N = 10;                                  % number of thruster firings
delta_T = T/N;                          % time between successive pulses

m = 0:N-1;
t_m = m*delta_T;                        % intervals when thrusters are fired
t_steps = length(t_m);

rho_m = zeros(1,t_steps);           % store distance to go values
rho_vecs = zeros(3,t_steps);        % store rho vectors (rho_m*rho_unit_vec)
rm = zeros(3,t_steps);              % change in distance to go from target
for i = 1:t_steps
    rho_m(i) = (rho_0*exp(a*t_m(i)) + (rho_dot_T/a)*(exp(a*t_m(i))-1)); % rho vector magnitude
    rho_vecs(:,i) = rho_unit_vec;           % populate with unit vectors
    rho_vecs(:,i) = rho_m(i)*rho_vecs(:,i); % multiply rho magnitude with rho unit vectors
    rm(:,i) = rT + rho_m(i)*rho_unit_vec;   % rm == rT at time T
end



plot3([0 r0(1)], [0 r0(2)], [0 r0(3)])  % r0 vector
hold("on")
plot3([0 rT(1)], [0 rT(2)], [0 rT(3)])  % rT vector
grid("on")

% 2d plot
% plot([0 r0(1)], [0 r0(3)])  % r0 vector
% hold("on")
% plot([0 rT(1)], [0 rT(3)])  % rT vector
% hold("on")
% plot([rT(1) r0(1)], [rT(3) r0(3)], "--")  % rho_0 vector
% hold("on")

% plot distance vectors
for i = 1:t_steps

    % 3d plot
    plot3([rT(1) rho_vecs(1,i)+rT(1)], [rT(2) rho_vecs(2,i)+rT(2)], [rT(3) ...
        rho_vecs(3,i)+rT(3)], "-", 'LineWidth',3)
    hold("on")

    % 2d plot
    % plot rho_vecs
%     plot([rT(1) rho_vecs(1,i)+rT(1)], [rT(3) rho_vecs(3,i)+rT(3)], "-", 'LineWidth',3)
%     hold("on")
    % plot rm
%     plot([0 rm(1,i)], [0 rm(3,i)], 'LineWidth',3)
%     hold("on")
end

grid("on")
xlabel("x position")
ylabel("y position")
zlabel("z position")


% parameters
mu = 398600;            % [km3/s2]
rRef = 6600;            % 222km altitude [km]
n = sqrt(mu/rRef^3);

t = delta_T;              % total trajectory time

% compute 
cw = [4-3*cos(n*t) 0 0 sin(n*t)/n 2*(1-cos(n*t))/n 0;
      6*(sin(n*t)-n*t) 1 0 2*(cos(n*t)-1)/n (4*sin(n*t)-3*n*t)/n 0;
      0 0 cos(n*t) 0 0 sin(n*t)/n;
      3*n*sin(n*t) 0 0 cos(n*t) 2*sin(n*t) 0;
      6*(n*cos(n*t)-n) 0 0 -2*sin(n*t) 4*cos(n*t)-3 0;
      0 0 -n*sin(n*t) 0 0 cos(n*t)];


% pull 3x3 matrices out of cw
A = cw(1:3,1:3);
B = cw(1:3,4:6);
C = cw(4:6,1:3);
D = cw(4:6,4:6);

% % compute the missing boundary conditions
% r = rm(:,2);        % final position from the target spacecraft
% v0 = B\(r-A*r0);    % compute the initial velocity to get to the final position
% v = C*r0 + D*v0;    % compute the final velocity at the final position
% 
% % sim
% deltaT = 0.01;
% x0 = [r0; v0];
% x_hist = cw_propagator(x0,n,t,deltaT);


% position vector of chaser relative to target from r0 to rT
r = rm;             % positions to arrive at
r(:,end+1) = rT;    % populate with the final arrival position

v0 = zeros(3,t_steps);  % storage for velocities to get to next waypoint
v = zeros(3,t_steps); % storage for velocities at the next waypoint

sim_delta_T = 0.01;     % resolution of CW sim
for i = 1:t_steps
    v0(:,i) = B\(r(:,i+1)-A*r(:,i));    % compute the initial velocity to get to the next waypoint
    v(:,i) = C*r(:,i) + D*v0(:,i);    % compute the final velocity at the next waypoint

    % sim
    x0 = [r(:,i); v0(:,i)];
    x_hist = cw_propagator(x0,n,t,sim_delta_T);

    % plotting trajectory
    x = x_hist(1,:);
    y = x_hist(2,:);
    z = x_hist(3,:);
    
    plot3(x,y,z)
    title('Glideslope Trajectory')
    xlabel('x position [km]') 
    ylabel('y position [km]')
    zlabel('z position [km]') 
    set(gca,'Ydir','reverse','Zdir','reverse')
    grid("on")
end

% plot velocity vector
% plot3([r0(1) v0(1,i)+r0(1)], [r0(2) v0(2,i)+r0(2)], [r0(3) ...
%     v0(3,i)+r0(3)], "-", 'LineWidth',3)
% hold("on")

% matches the magnitude of rho_dot_0 = -1, rho_dot_T = -0.3
% i.e. generated velocity vectors matches the desired inital and final
% constraints
rho_dot_T = -0.3;
v0_norm = norm(v(:,1))
vT_norm = norm(v(:,end))
