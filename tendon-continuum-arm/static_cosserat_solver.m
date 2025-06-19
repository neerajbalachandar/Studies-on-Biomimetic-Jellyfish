(* Static Cosserat Rod Solver *)
(* Constant Kse,Kbt, Converging Tendon Routing *)
(* Varying Tendon Tensions, Neglecting external static forces *)

__author__ = "Neeraj Balachandar"
__date__ = "19 June 2025"
__contact__ = "neerajbalachandar@gmail.com"


params = struct();

params.h_initial = 0.005;
params.L_total = 0.8; 

% Material properties (stiffness matrices)
E = 79e6;                % Young's modulus [Pa]
G_shear = 27e6;                 % Shear modulus [Pa]
A = pi*(0.001)^2;         % Cross-sectional area [m²]
I = pi*(0.001)^4/4;       % Second moment of area [m⁴]
J = 2*I;                     % Torsional constant [m⁴]
params.Kse = diag([G_shear*A, G_shear*A, E*A]); % Shear/extension stiffness
params.Kbt = diag([E*I, E*I, G_shear*J]); % Bending/torsion stiffness

% Reference configuration
params.v_star = [0; 0; 1];          % Undeformed shear/stretch
params.u_star = [0; 0; 0];          % Undeformed curvature/twist

% External loads (gravity example)
rho = 1200;                         % Density [kg/m³]
g = 9.81;                           % Gravity [m/s²]

% Initial state (straight rod)
p0 = [0; 0; 0];                     % Initial position [m]
R0 = [0, 0, 1;                      
      0, 1, 0;
     1, 0, 0];                     
v0 = params.v_star;                 % Initial strain
u0 = params.u_star;                 % Initial curvature
y0 = [p0; R0(:); v0; u0];           % Pack initial state

% Solve ODE from base (s=0) to tip (s=0.5m)
s_span = [0, 1];                  % Arc length range [m]

tau_values = 0:0.05:0.2;
solutions = cell(length(tau_values), 1);


% Solve for each tension value
for i = 1:length(tau_values)
    params.tau = tau_values(i);
    solutions{i} = cosserat_static_solver(s_span, y0, params);
end

% Plot all configurations
plot_multiple_rod_shapes(solutions, tau_values);


%--------------------FUNCTION DEFINITIONS------------------------
function [r, r_dot] = compute_tendon_routing(s, params)
    
    h = params.h_initial;           
    L = params.L_total;            

    % Based on the relationship: y = h - h*s/L
    convergence_factor = 1 - s/L;   % Goes from 1 to 0 as s goes from 0 to L
    
    % Routing vector r(s) - tendon position in body frame
    r = [convergence_factor * h;    % x-component: converges to 0
         0;                         % y-component: stays at 0
         0];                        % z-component: along centerline
    
    r_dot = [-h/L;            
             0;             
             0];               
end



function y_out = cosserat_static_solver(s_span, y0, params)
    % Solves the Cosserat rod ODEs statically for a tendon-driven continuum robot
    % s_span: [s_start, s_end]
    % y0: initial state vector [p0; R0(:); v0; u0]
    % params: struct with robot and tendon parameters

    % ODE integration
    [s, y] = ode45(@(s, y) cosserat_rod_rhs(s, y, params), s_span, y0);

    y_out.s = s;
    y_out.y = y;
end

function dy_ds = cosserat_rod_rhs(s, y, params)
    % Unpack state
    p = y(1:3);                    % position (3x1)
    R = reshape(y(4:12), 3, 3);    % rotation matrix (3x3)
    v = y(13:15);                  % strain vector (3x1)
    u = y(16:18);                  % curvature vector (3x1)

    % Compute system matrices and RHS (user-defined functions)
    [Kse, Kbt, A, B, G, H, d, c] = compute_system_matrices(s, p, R, v, u, params);

    % Assemble block matrix and invert
    M = [Kse + A, G; B, Kbt + H];
    rhs = [d; c];
    strain_dot = M \ rhs;   % Equivalent to inv(M)*rhs, but more efficient

    % Kinematics
    p_dot = R * v;
    u_hat = hat(u);         % Skew-symmetric matrix from u
    R_dot = R * u_hat;

    % Output derivative of state
    dy_ds = zeros(18, 1);
    dy_ds(1:3) = p_dot(:);
    dy_ds(4:12) = R_dot(:);
    dy_ds(13:18) = strain_dot;
end

function [Kse, Kbt, A, B, G, H, d, c] = compute_system_matrices(s, p, R, v, u, params)
    % Compute converging tendon routing at current arc length
    [r, r_dot] = compute_tendon_routing(s, params);
    
    % Parameters for converging tendon
    T = params.tau;                 % Tendon tension (scalar)
    Kse = params.Kse;              % Shear/extension stiffness (3x3)
    Kbt = params.Kbt;              % Bending/torsion stiffness (3x3)
    v_star = params.v_star;        % Reference strain (3x1)
    u_star = params.u_star;        % Reference curvature (3x1)
    
    % Hydrostatic forces (from previous discussion)
    %ho_water = 1025;              % Seawater density [kg/m³]
    %ho_robot = 8000;              % Robot density [kg/m³]
    g = 9.81;                      % Gravity [m/s²]
    A_cross = pi*(0.001)^2;        % Cross-sectional area [m²]

    f_external = [0;0;0];
    l_external = [0;0;0];
    
    %_external_global = (rho_water - rho_robot) * g * A_cross * [0; 0; 0.8];
    %_external_local = R' * f_external_global;  
    %
    %r_local = [0;0;0];
    %_external_local = cross(r_local, f_external_local);
    %_external_global = R * l_external_local;

    % Skew-symmetric matrices
    u_hat = hat(u);
    r_hat = hat(r);
    
    % Tendon path derivative in body frame
    p_dot_b = u_hat * r + r_dot + v;
    norm_p_dot_b = norm(p_dot_b);
    
    % Avoid division by zero
    if norm_p_dot_b < 1e-8
        p_dot_b_unit = [0; 0; 1];
        norm_p_dot_b = 1;
    else
        p_dot_b_unit = p_dot_b / norm_p_dot_b;
    end
    
    % System matrices for converging tendon
    A = -T * (p_dot_b_unit * p_dot_b_unit');
    B = r_hat * A;
    G = -A * r_hat;
    H = -B * r_hat;
    
    % RHS vectors
    a = A * (u_hat * p_dot_b + u_hat * r + r_dot);
    b = r_hat * a;
    
    % Full RHS vectors
    d = Kse * v_star - u_hat * Kse * (v - v_star) - f_external - a;
    c = Kbt * u_star - u_hat * Kbt * (u - u_star) -  l_external - b;
end


function S = hat(v)
    S = [  0   -v(3)  v(2);
          v(3)   0   -v(1);
         -v(2) v(1)    0  ];
end



%-------------------------PLOTTING-------------------------------
function plot_multiple_rod_shapes(solutions, tau_values)
    figure;
    hold on;
    % Use viridis colormap (perceptually uniform, colorblind-friendly)
    colors = viridis(length(tau_values));
    % Plot configurations with depth cues
    for i = 1:length(solutions)
        pos = solutions{i}.y(:,1:3);
        plot3(pos(:,1), pos(:,2), pos(:,3), ...
              'Color', colors(i,:), 'LineWidth', 2.5, ...
              'LineStyle', '-', 'Marker', 'none');
    end
    
    % Enhanced marker styling
    pos = solutions{end}.y(:,1:3);
    scatter3(pos(1,1), pos(1,2), pos(1,3), 40, ...
            'MarkerFaceColor', [0.8 0.4 0.4], ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.5);
    
    scatter3(pos(end,1), pos(end,2), pos(end,3), 40, ...
            'MarkerFaceColor', [0.4 0.4 0.8], ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.5); 

    ax = gca;
    ax.FontName = 'LaTeX';
    ax.FontSize = 12;
    ax.XLabel.String = 'X [m]';
    ax.YLabel.String = 'Y [m]';
    ax.ZLabel.String = 'Z [m]';
    ax.GridColor = [0.4 0.4 0.4];
    ax.GridAlpha = 0.3;
    ax.XColor = [0.3 0.3 0.3];
    ax.YColor = [0.3 0.3 0.3];
    ax.ZColor = [0.3 0.3 0.3];
    
    title('Static Model of Tendon Driven Continuum Lappet using Coupled Cosserat Rod and String Theory for Varying Tendon Converging Configuration Tension', ...
          'FontSize', 14, 'FontWeight', 'bold');

    legend(arrayfun(@(t) sprintf('$\\tau = %d$ N', t), tau_values, ...
    'UniformOutput', false), ...
    'Interpreter', 'latex', ...
    'Location', 'northeastoutside');

    % View and lighting
    view(3);
    axis equal tight;
    grid on;
    box on;
    light('Style','infinite','Position',[1 0.5 1]);
    lighting gouraud
    material dull
    hold off;
end
% Helper function for viridis colormap (matches search results[1])
function colors = viridis(n)
    % Viridis RGB values from search result[1]
    viridis_map = [0.267004, 0.004874, 0.329415
                   0.267968, 0.223549, 0.512008
                   0.190631, 0.407061, 0.556089
                   0.127568, 0.566949, 0.550556
                   0.20803, 0.718701, 0.472873
                   0.565498, 0.84243, 0.262877
                   0.993248, 0.906157, 0.143936];
    
    % Interpolate for smooth color transitions
    x = linspace(0, 1, size(viridis_map,1));
    xi = linspace(0, 1, n);
    colors = [interp1(x, viridis_map(:,1), xi)', ...
              interp1(x, viridis_map(:,2), xi)', ...
              interp1(x, viridis_map(:,3), xi)'];
end

