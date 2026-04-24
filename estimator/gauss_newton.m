function [R_opt, bias_opt] = gauss_newton(data, lambda, R_init)
% Estimate complete 3D rotation matrix and clock drift bias using Gauss-Newton iteration based on Doppler measurements
%
% Input:
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
% Output:
%   R_opt: optimized rotation matrix
%   bias_opt: clock drift (m/s)

max_iter = 20;
eps = 1e-6;

R = R_init;
bias = 0;

for iter = 1:max_iter
    J_all = [];
    r_all = [];
    
    for k = 1:length(data)
        obs = data{k};
        local_vel = obs.local_velocity; % 3x1 local velocity
        
        N = length(obs.visible_sats);
        for i = 1:N
            s = obs.visible_sats(i);
            sat_pos = obs.sat_positions(:, s);
            sat_vel = obs.sat_velocities(:, s);
            rx_pos = obs.ecef_position;
            doppler_meas = obs.doppler(s);
            
            los = sat_pos - rx_pos;
            los = los / norm(los);
            
            % Predicted Doppler (unit m/s, note that lambda is already multiplied here)
            doppler_pred = los' * (sat_vel - R * local_vel) + bias;
            res = doppler_meas - doppler_pred;

            % Cross product matrix
            skew_local_vel = hat(R * local_vel); % hat operation
              
            J_rot =- los' * skew_local_vel; % 1x3
      
            % Partial derivative with respect to bias
            J_bias = 1;
            
            J_all = [J_all; J_rot, J_bias];
            r_all = [r_all; res];
        end
    end
    
    % Gauss-Newton update
    dx = (J_all' * J_all) \ (J_all' * r_all);
    delta_w = dx(1:3);
    delta_bias = dx(4);
    
    % Rotation matrix update
    R = expm(hat(delta_w))'*R ;
    bias = bias + delta_bias;
    
    % Hess = J_all' * J_all;
    % [U,S,~] = svd(Hess);
    % s = diag(S)

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
