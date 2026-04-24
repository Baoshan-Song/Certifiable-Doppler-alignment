%% Parameter setup
rng(4);
motion_type = 'step3D';  % 'step3D' or 'spiral'

dt = 1;                % Time interval [s]
N =10;                % Simulation duration [step]
c = 299792458;         % Speed of light [m/s]
f_L1 = 1.57542e9;      % L1 frequency [Hz]
lambda = c / f_L1;     % Wavelength [m]
ecef_start = [3980581; 24; 4966824];  % Initial ECEF position
mean_velocity =3;
% Satellite parameters
M =4;      % Number of orbital planes
P =4;      % Satellites per plane
num_sats = M * P;
i_deg = 55;
i_rad = deg2rad(i_deg);
r_orbit = 2.656e7;
omega = 2*pi / (12*3600);

%% Generate trajectory
v_ecef_true_all = zeros(3, N);
ecef_positions = zeros(3, N);

for i = 1:N
    switch motion_type
        case 'step3D'
            if i <= N/3*2
                v = [mean_velocity; 0; 0];
            else
                v = [0; 0; mean_velocity];
            end
        case 'spiral'
            theta = 2*pi*(i-1)/N;
            v = mean_velocity * [cos(theta); sin(theta);0];
        otherwise
            error('Unknown motion type');
    end

    v_ecef_true_all(:, i) = v;
    if i == 1
        ecef_positions(:, i) = ecef_start;
    else
        ecef_positions(:, i) = ecef_positions(:, i-1) + dt * v;
    end
end


% Walker Delta constellation + random orbital phase
sat_positions = zeros(3, num_sats, N);
sat_velocities = zeros(3, num_sats, N);

for m = 1:M
    for p = 1:P
        s = (m-1)*P + p;

        % RAAN for each orbital plane
        RAAN = 2*pi * (m-1)/M;

        % Initial position of each satellite on the orbit is random
        initial_phase = 2*pi * rand();

        % Construct orbital rotation matrix
        R_inc = [1 0 0; 0 cos(i_rad) -sin(i_rad); 0 sin(i_rad) cos(i_rad)];
        R_raan = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
        R_total = R_raan * R_inc;

        for k = 1:N
            t = (k-1)*dt;
            theta = omega * t + initial_phase;

            % Position and velocity in orbital plane
            pos_orb = r_orbit * [cos(theta); sin(theta); 0];
            v_orb = r_orbit * omega * [-sin(theta); cos(theta); 0];

            % Transform to Earth-centered inertial coordinate system
            sat_positions(:, s, k) = R_total * pos_orb;
            sat_velocities(:, s, k) = R_total * v_orb;
        end
    end
end


%% Receiver orientation (R_ecef2local)
R_ecef2local = random_rot_matrix(3);

%% Doppler measurement calculation
data_cells = cell(1, N);
doppler_measurements = zeros(num_sats, N);
rev_clk_drift = 1000;

for i = 1:N
    rx_pos = ecef_positions(:, i);
    rx_vel = v_ecef_true_all(:, i);

    doppler = nan(num_sats, 1);
    sat_pos_list = nan(3, num_sats);
    sat_vel_list = nan(3, num_sats);

    for s = 1:num_sats
        sat_pos = sat_positions(:, s, i);
        sat_vel = sat_velocities(:, s, i) ;

        los_ecef = sat_pos - rx_pos+100*rand(1);
        los_unit = los_ecef / norm(los_ecef);

        los_local = R_ecef2local * los_unit;
        el = asin(los_local(3));
        if rad2deg(el) < 10
            continue;
        end

        rel_vel = sat_vel - rx_vel;% + rand(1);%*1;
        doppler_true = dot(rel_vel, los_unit)/lambda + rev_clk_drift+ rand(1)*5;
        
        doppler(s) = doppler_true * lambda;  % Convert back to [m/s]

        sat_pos_list(:, s) = sat_pos;
        sat_vel_list(:, s) = sat_vel;
    end

    data.time = (i-1)*dt;
    data.ecef_position = rx_pos;
    data.local_velocity = R_ecef2local * rx_vel;
    data.sat_positions = sat_pos_list;
    data.sat_velocities = sat_vel_list;
    data.doppler = doppler;
    data.true_attitude = R_ecef2local;
    data.visible_sats = find(~isnan(doppler));

    data_cells{i} = data;
    doppler_measurements(:, i) = doppler;
end

%% Visualization
figure;
plot3(ecef_positions(1,:), ecef_positions(2,:), ecef_positions(3,:), 'b.-');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Receiver ECEF Trajectory');
grid on; axis equal;

figure;
plot((0:N-1)*dt, doppler_measurements');
xlabel('Time [s]'); ylabel('Doppler [m/s]');
title('Doppler Measurements');
legend(arrayfun(@(x) sprintf('Sat %d', x), 1:num_sats, 'UniformOutput', false));
grid on;

