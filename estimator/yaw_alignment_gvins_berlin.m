function [yaw, bias, R] = yaw_alignment_gvins_berlin(data, lambda, R_e_n, init_yaw)
% YAW_ALIGNMENT_GVINS Estimate yaw and clock bias
% 
% Inputs:
%   data: struct for each frame, contains following fields
%       .sat_positions: 3×N satellite positions
%       .sat_velocities: 3×N satellite velocities
%       .ecef_position:  3×1 receiver position
%       .doppler: 1×N Doppler measurements (Hz)
%   lambda: wavelength (unit: m)
%   R_rp: known rotation matrix for roll and pitch (3x3), ECEF → body_rp
%
% Outputs:
%   yaw: horizontal angle in rad
%   bias: clock bias (unit: m/s)
%   R: total rotation R = R_rp * Rz(yaw), i.e., ECEF → body

    % Initial value
    x = [init_yaw; 0]; % [yaw; bias]

    max_iter = 2000;
    eps = 1e-6;

    for iter = 1:max_iter
        J = [];
        r = [];

        yaw = x(1);
        bias = x(2);

        % test --sbs
        pitch = -0 ;roll  =-0 ;
        Rz = [cos(yaw), -sin(yaw), 0;
              sin(yaw),  cos(yaw), 0;
              0,         0,        1];

    % Rz = [cos(yaw) -sin(yaw) 0;
    %       sin(yaw)  cos(yaw) 0;
    %       0         0        1];
    Ry = [cos(pitch) 0 sin(pitch);
          0          1 0;
         -sin(pitch) 0 cos(pitch)];
    Rx = [1 0 0;
          0 cos(roll) -sin(roll);
          0 sin(roll)  cos(roll)];
    % 
    % R_n_b = Rz * Ry * Rx;
        R_n_b =  Rz* Ry * Rx;

        for k = 1:size(data.gps_week, 1)
            obs = data.sat_position{k};
            local_vel =  data.local_velocity(:, k);
            % N = size(obs.sat_positions, 2);
            N =size(obs, 2);
            for i = 1:N
                % sat_idx = obs.visible_sats(i);
                los = data.sat_position{k}(:, i)- data.ecef_position(k, :)';
                los = los / norm(los);
                vs = data.sat_velocity{k}(:, i);
                doppler =-data.doppler{k}(i);

                % Predicted velocity projection (unit: m/s)
                % pred = los' * R_total * v_s;
                % resi = -doppler + pred + bias;

                doppler_predicted =los'* (vs - R_e_n*R_n_b*local_vel) + bias;
                resi = doppler - doppler_predicted;

                % Construct Jacobian
                dRz_dyaw = [-sin(yaw), -cos(yaw), 0;
                             cos(yaw), -sin(yaw), 0;
                             0,         0,        0];
                dR_dyaw = R_e_n* dRz_dyaw* Ry * Rx;

                % dpred_dyaw = los' * dR_dyaw * v_s;
                % dpred_dbias = 1;

                J_yaw = -los'* dR_dyaw*local_vel;
                J_bias = 1;


                % J = [J; dpred_dyaw, dpred_dbias];
                % r = [r; resi];
                J = [J; J_yaw, J_bias];
                r = [r; resi];

            end
        end

        % Gauss-Newton step
        dx = (J' * J) \ (J' * r);
        x = x + dx;

        if norm(dx) < eps
            break;
        end
    end


    aligned_yaw = mod(x(1)+pi, 2*pi) - pi;

    yaw = aligned_yaw;
    bias = x(2);
    R =  R_e_n* [cos(yaw), -sin(yaw), 0;
                sin(yaw),  cos(yaw), 0;
                0,         0,        1] * Ry * Rx;
end
