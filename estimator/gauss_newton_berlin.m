function [R_opt, bias_opt] = gauss_newton_berlin(data, lambda, R_init)
% Estimate complete 3D rotation matrix and clock bias using Gauss-Newton iteration based on Doppler measurements
%
% Inputs:
%   data: cell array, each frame contains
%       .sat_positions: 3×N satellite positions (ECEF)
%       .sat_velocities: 3×N satellite velocities (ECEF)
%       .ecef_position: 3×1 receiver position (ECEF)
%       .doppler: 1×N Doppler measurements (Hz)
%       .local_velocity: 3×1 local velocity (local frame)
%       .visible_sats: visible satellite indices
%   lambda: wavelength, unit m
%   R_init: initial rotation matrix, ECEF→local
%
% Outputs:
%   R_opt: optimized rotation matrix
%   bias_opt: clock bias (m/s)

max_iter = 10;
eps = 1e-6;

R = R_init;
bias = 0;

down_sample_cnt =0;


for iter = 1:max_iter
    J_all = [];
    r_all = [];
    W = [];
    for k = 1:size(data.gps_week, 1)

        % 
        if mod(k,10)~=0
            continue;
        end


        obs = data.sat_position{k};
        local_vel = data.local_velocity(:, k); % 3x1 local velocity

           sat_thres =4;
        N =size(obs, 2); 
      max_sat = 0;
      if(N<sat_thres)
          max_sat = N;
      else
          max_sat =sat_thres;
      end

        for i = 1:max_sat % --sbs
            sat_pos = data.sat_position{k}(:, i);
            sat_vel = data.sat_velocity{k}(:, i);
            rx_pos = data.ecef_position(k, :)';
            doppler_meas =-data.doppler{k}(i);
            
            % --sbs
            if isnan(doppler_meas)
                continue;
            end

            los = sat_pos - rx_pos;
            los = los / norm(los);
            
            % Predicted Doppler (unit: m/s, note that lambda is already multiplied here)
            doppler_pred = los' * (sat_vel - R * local_vel) + bias;
            % doppler_pred = los' * (sat_vel - R * local_vel) + bias;
            res = doppler_meas - doppler_pred;
            
            % Simple elevation angle dependent weight model
            varD90 = 0.1^2; % (m/s)^2
            el = satElevationAngle(sat_pos, rx_pos);
            w = 1./(varD90./sin(el)); 
            W = [W;w];

            % Cross product matrix
            skew_local_vel = hat(R * local_vel); % hat operation
  
            J_rot =- los' * skew_local_vel; % 1x3

            J_bias = 1;
            
            J_all = [J_all; J_rot, J_bias];
            r_all = [r_all; res];
        end
    end
    W = diag(W);
    % Gauss-Newton update
    dx = (J_all' * J_all) \ (J_all' * r_all);
    delta_w = dx(1:3);
    delta_bias = dx(4);
    
    % Rotation matrix update
    R = expm(hat(delta_w))'*R ;
    bias = bias + delta_bias;
    
    H = J_all' * J_all;
    eigvals = eig(H);
    min_eig = min(eigvals);

    if norm(dx) < eps
        H = J_all' * J_all;
        eigvals = eig(H);
        min_eig = min(eigvals);
        % disp(['H minimum eigenvalue is: ', num2str(min_eig)]);
        break;
    end
end

R_opt = R;
bias_opt = bias;

end

function S = hat(v)
S = [  0   -v(3)  v(2);
      v(3)  0   -v(1);
     -v(2) v(1)   0];
end



function el = satElevationAngle(XYZ_sat, XYZ_rcv)
% satElevationAngle Compute satellite elevation angle relative to receiver
% Inputs:
%   XYZ_sat - satellite ECEF coordinates [3x1] (m)
%   XYZ_rcv - receiver ECEF coordinates [3x1] (m)
% Outputs:
%   el - elevation angle (rad)

    % WGS84 parameters
    a = 6378137.0;          % semi-major axis
    f = 1/298.257223563;    % flattening
    e2 = f*(2-f);           % square of first eccentricity

    % Receiver ECEF -> BLH
    x = XYZ_rcv(1); y = XYZ_rcv(2); z = XYZ_rcv(3);
    lon = atan2(y,x);
    p = sqrt(x^2 + y^2);
    lat = atan2(z, p*(1-e2)); % initial value
    for k = 1:5  % Newton iteration to correct latitude
        N = a / sqrt(1 - e2*sin(lat)^2);
        h = p/cos(lat) - N;
        lat = atan2(z + e2*N*sin(lat), p);
    end

    % ECEF -> ENU transformation matrix
    R = [ -sin(lon),            cos(lon),           0; ...
          -sin(lat)*cos(lon),  -sin(lat)*sin(lon),  cos(lat); ...
           cos(lat)*cos(lon),   cos(lat)*sin(lon),  sin(lat)];

    % Satellite vector (ECEF)
    rho = XYZ_sat - XYZ_rcv;

    % Convert to ENU
    enu = R * rho;

    % Calculate elevation angle
    east = enu(1); north = enu(2); up = enu(3);
    el = atan2(up, sqrt(east^2 + north^2));
end
