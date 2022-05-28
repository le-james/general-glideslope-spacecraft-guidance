%% mission design
clear; clc;

number_of_targets = 3;
mission_parameters = 7;
% includes initial chaser position
mission_targets = zeros(mission_parameters,number_of_targets+1);

% choosen by mission requirements
r0 = [500 0 -200]'; % chaser inital position relative to the target [m] - 3x1 vector

rho_dot_0 = -1;    % inital commanded velocity [m/s]
rho_dot_T = -0.3;  % velocity to arrive at target position [m/s]

%{
    mission targets (user defined parameters) - 7x1 vector
    mt_0(1:3)==r0: initial cartesian position of chaser relative to target
                   position mt_#+1(1:3) - (3x1 vector)
    mt_T#(1:3)==rT: goal cartesian position to move the chaser to -
                     (3x1 vector)
    mt_Tn: the last target postion - mt_Tn(4:7) are left as zeros
    mt_0(4) and mt_T#(4): magnitude of initial desired delta V to move to 
                          the next target position mt_#+1(1:3)
    mt_0(5) and mt_T#(5): magnitude of desired velocity to arrive at 
                          mt_#+1(1:3) from mt_#(1:3)
    mt_0(6) and mt_T#(6): desired number of thruster firing firings at 
                          mt_#(1:3) to get to mt_#+1(1:3)
    mt_0(7) and mt_T#(7): maneuver type: 0 = inbound, 1 = outbound, 
                          2 = circumnav from mt_#(1:3) to mt_#+1(1:3)
%}
mt_0 = [r0(1) r0(2) r0(3) rho_dot_0 rho_dot_T 10 0]';  % initial
mt_T1 = [200 0 -200 rho_dot_0 rho_dot_T 10 0]';   
mt_T2 = [75 0 -100 rho_dot_0 rho_dot_T 10 0]';
mt_Tn = [50 0 0 0 0 0 0]';

% store designed mission - order of chaser goals positions - 7x1 vector
mission_targets(:,1) = mt_0;
mission_targets(:,2) = mt_T1;
mission_targets(:,3) = mt_T2;
mission_targets(:,4) = mt_Tn;

%% generate glideslope guidance trajectory

% storage of parameters for trajectory generation - described below and in
% reference paper
rho_0 = zeros(1,number_of_targets);
cos_alpha = zeros(1,number_of_targets);
cos_beta = zeros(1,number_of_targets);
cos_gamma = zeros(1,number_of_targets);
rho_unit_vec = zeros(3,number_of_targets);
a = zeros(1,number_of_targets);
T = zeros(1,number_of_targets);
delta_T = zeros(1,number_of_targets);
m = zeros(max(mission_targets(6,:)),number_of_targets);
t_m = zeros(max(mission_targets(6,:)),number_of_targets);

rho_m = zeros(max(mission_targets(6,:)),number_of_targets);      % store distance to go values
rho_vecs = zeros(3,max(mission_targets(6,:)),number_of_targets); % store rho vectors (rho_m*rho_unit_vec)
rm = zeros(3,max(mission_targets(6,:)),number_of_targets);       % change in distance to go from target

for i = 1:number_of_targets

    r0 = mission_targets(1:3,i);    % initial chaser position
    rT = mission_targets(1:3,i+1);  % chaser target position

    % magnitude of rho_0
    rho_0(i) = norm(r0-rT);

    % pull out inital and target positions
    x0 = r0(1);
    y0 = r0(2);
    z0 = r0(3);
    xT = rT(1);
    yT = rT(2);
    zT = rT(3);

    % compute direction cosines to compute the unit rho unit vector
    cos_alpha(i) = (x0-xT)/rho_0(i);
    cos_beta(i) = (y0-yT)/rho_0(i);
    cos_gamma(i) = (z0-zT)/rho_0(i);

    % rho unit direction vector
    rho_unit_vec(1:3,i) = [cos_alpha(i) cos_beta(i) cos_gamma(i)]';
    
    % a<0
    a(i) = (mission_targets(4,i)-mission_targets(5,i))/rho_0(i);

    % total time of rendezvous trajectory
    T(i) = (1/a(i))*log(mission_targets(5,i)/mission_targets(4,i));

    % time between successive pulses
    delta_T(i) = T(i)/mission_targets(6,i);

    % number of thruster firings (pulses) to move to the goal position
    m(1:mission_targets(6,i),i) = (0:mission_targets(6,i)-1)';

    % intervals when thrusters are fired
    t_m(:,i) = m(:,i)*delta_T(i);

    % need to figure this part out
    for j = 1:mission_targets(6,i)
        
        % rho vector magnitude
        rho_m(j,i) = (rho_0(i)*exp(a(i)*t_m(j,i)) + (mission_targets(5,i)/a(i))*(exp(a(i)*t_m(j,i))-1));
        
        % populate with unit vectors
        rho_vecs(:,j,i) = rho_unit_vec(1:3,i);
        
        % multiply rho magnitude with rho unit vectors
        rho_vecs(:,j,i) = rho_m(j,i)*rho_vecs(:,j,i);
        
        % rm == rT at time T
        rm(:,j,i) = mission_targets(1:3,i) + rho_m(j,i)*rho_unit_vec(1:3,i);
    end
end

%% for plots
% labels
labels = strings(1,number_of_targets+2);
labels(1) = 'Chaser S/C';
labels(end) = 'Target S/C';

% store chaser initial positions, target positions and target s/c positiion
points = mission_targets(1:3,:);
r0 = points(:,1);
rn = points(:,end);
% for plots

%% figure 1 plot of mission goals key points
figure
% plot initial position of the chaser
plot3(r0(1),r0(2),r0(3),"o",'LineWidth',2);hold("on")  % r0 vector

% plot mission targets
for i = 1:number_of_targets

    % label target goal points
    labels(i+1) = sprintf('Target %d',i);

    plot3(points(1,i+1),points(2,i+1),points(3,i+1),"o",'LineWidth',2);hold("on")
end

% % add target spacecraft label to plot
% text(0,0,0,'Target S/C','VerticalAlignment','bottom','HorizontalAlignment','right')
% % add chaser initial position label to plot
% text(r0(1),r0(2),r0(3),'Initial Chaser S/C','VerticalAlignment','bottom','HorizontalAlignment','left')
% % add target goal points label to plot
% text(mission_targets(1,:),mission_targets(2,:),mission_targets(3,:),labels, ...
%     'VerticalAlignment','bottom','HorizontalAlignment','left')

plot3(0,0,0,"o",'LineWidth',2);hold("on")  % target s/c

% add labels to plot
legend(labels)

title('plot')
xlabel('x Position [m]')
ylabel('y Position [m]')
zlabel('z Position [m]')
grid("on")

%% figure 2 plot of glideslope from chaser to target goal points
figure
% plot initial position of the chaser
plot3(r0(1),r0(2),r0(3),"o",'LineWidth',2);hold("on")  % r0 vector

% plot the glideslope
for i = 1:number_of_targets

    % label target goal points
    labels(i+1) = sprintf('Target %d',i);
    
%     quiver3(points(1,i),points(2,i),points(3,i), ...
%     points(1,i+1)-points(1,i),points(2,i+1)-points(2,i),points(3,i+1)-points(3,i), ...
%     'off','-.','LineWidth',2,'MaxHeadSize',999)
    
    for j = 1:length(mission_targets(6,i))
        plot3([points(1,i+1) points(1,i+1)+rho_vecs(1,j,i)], ...
            [points(2,i+1) points(2,i+1)+rho_vecs(2,j,i)], ...
            [points(3,i+1) points(3,i+1)+rho_vecs(3,j,i)], ...
            '-','LineWidth',2); hold("on")
%     quiver3(points(1,i),points(2,i),points(3,i), ...
%         points(1,i+1)-points(1,i),points(2,i+1)-points(2,i),points(3,i+1)-points(3,i), ...
%         'off','-.','LineWidth',2,'MaxHeadSize',999)
    end



%     % plot the glideslope from
%     plot3([points(1,i) points(1,i+1)],[points(2,i) points(2,i+1)], ...
%         [points(3,i) points(3,i+1)],'--','LineWidth',2)

%     % glideslope arrow plots
%     quiver3(points(1,i),points(2,i),points(3,i), ...
%         points(1,i+1)-points(1,i),points(2,i+1)-points(2,i),points(3,i+1)-points(3,i), ...
%         'off','-.','LineWidth',2,'MaxHeadSize',999)

end

% didn't finish dev - ploting the rho vector instead of using this
% for j = 1:length(points)-1
%     vec_sub = points(:,j+1) - points(:,j); % vector difference between points
%     vec_norm = norm(vec_sub);              % take 2-norm of vector 
%     vec_unit = vec_sub/vec_norm;           % convert to unit vectors
%     
%     % one vector to another spaced by the unit vectors in between them
%     line_field = [];
%     if vec_sub(1) == 0 && vec_sub(2) ~= 0 && vec_sub(3) ~= 0
%         line_field(2,:) = points(2,j):vec_unit(2):points(2,j+1);
%         line_field(3,:) = points(3,j):vec_unit(3):points(3,j+1);
%     elseif
%     elseif
%     else
%     end
%     line_field = [points(1,j):vec_unit(1):points(1,j+1);
%                   points(2,j):vec_unit(2):points(2,j+1);
%                   points(3,j):vec_unit(3):points(3,j+1)];
% 
%     quiver3(line_field(1,j),line_field(2,j),line_field(3,j), ...
%     line_field(1,j+1)-line_field(1,j),line_field(2,j+1)-line_field(2,j), ...
%     line_field(3,j+1)-line_field(3,j), ...
%     'off','-.','LineWidth',2,'MaxHeadSize',999)
%     hold("on")
% end

plot3(0,0,0,"o",'LineWidth',2);hold("on") % target s/c

% add labels to plot
legend(labels)

title('Glideslope Vectors')
xlabel('x Position [m]')
ylabel('y Position [m]')
zlabel('z Position [m]')
grid("on")

