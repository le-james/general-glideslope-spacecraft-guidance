function x_hist = cw_propagator(x0,n,t,deltaT)

% Returns the rendezvous trajectory of a chaser and target spacecraft in 
% circular orbit.
%   x_hist = cw_propagator(x0,n,t,deltaT)
%   INPUTS:
%     r0: 6x1 position and velocity vector of the chaser w.r.t to the target [km]
%     n: mean motion of target [1/s]
%     t: time of rendezvous trajectory [s]
%     t_inc: simulation time steps, smaller number = finer resolution [s]
%   OUTPUTS:
%     x_hist: trajectory of rendezvous [km]

%     v: velocity at final position (dock with target), fire v in negative 
%        direction to get relative velocity of zero [km/s]


    % simulation times
    t_hist = 0:deltaT:t;

    % total simulation steps
    t_steps = length(t_hist);


    % clohessy-whiltshire equations in matrix form
    function [xn,cw] = cw_equations(x0,n,t)
        cw = [4-3*cos(n*t) 0 0 sin(n*t)/n 2*(1-cos(n*t))/n 0;
              6*(sin(n*t)-n*t) 1 0 2*(cos(n*t)-1)/n (4*sin(n*t)-3*n*t)/n 0;
              0 0 cos(n*t) 0 0 sin(n*t)/n;
              3*n*sin(n*t) 0 0 cos(n*t) 2*sin(n*t) 0;
              6*(n*cos(n*t)-n) 0 0 -2*sin(n*t) 4*cos(n*t)-3 0;
              0 0 -n*sin(n*t) 0 0 cos(n*t)];

        xn = cw*x0;  % compute the next states
    end

%   USED FOR COMPUTING THE INVERSE PROBLEM OF CW - given r0 find v0
%     % cw matrix at time t
%     [~, cwBound] = cw_equations(zeros(6)',n,t);
% 
%     % pull 3x3 matrices out of cw
%     A = cwBound(1:3,1:3);
%     B = cwBound(1:3,4:6);
%     C = cwBound(4:6,1:3);
%     D = cwBound(4:6,4:6);
%     
%     % compute the missing boundary conditions
%     r = [0 0 0]';       % final position from the target spacecraft
%     v0 = B\(r-A*r0);    % compute the initial velocity to get to the final position
%     v = C*r0 + D*v0;    % compute the final velocity at the final position
%     
%     x0 = [r0; v0];      % initial states


    % simulatation
    x_hist = zeros(length(x0), t_steps);    % storage for time history of states
%     x_hist(:,1) = x0;                       % store initial states                       
    for i = 1:t_steps                                    
        [xn, ~] = cw_equations(x0, n, t_hist(i));  
        x_hist(:,i) = xn;
    end

%     % plotting trajectory
%     x = x_hist(1,:);
%     y = x_hist(2,:);
%     z = x_hist(3,:);
%     
%     plot3(x,y,z)
%     title('Glideslope Trajectory')
%     xlabel('x position [km]') 
%     ylabel('y position [km]')
%     zlabel('z position [km]') 
%     grid("on")
end