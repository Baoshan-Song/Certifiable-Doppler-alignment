function [yaw, bias, R] = yaw_alignment_gvins(data, lambda, R_rp, init_yaw)
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

    max_iter = 20;
    eps = 1e-6;

    for iter = 1:max_iter
        J = [];
        r = [];

        yaw = x(1);
        bias = x(2);
        
        Rz = [cos(yaw), -sin(yaw), 0;
              sin(yaw),  cos(yaw), 0;
              0,         0,        1];
        R_total = R_rp * Rz;

        for k = 1:length(data)
            obs = data{k};
            local_vel = obs.local_velocity;
            % N = size(obs.sat_positions, 2);
            N = size(obs.visible_sats);
            for i = 1:N
                sat_idx = obs.visible_sats(i);
                los = obs.sat_positions(:,sat_idx)-obs.ecef_position;
                los = los / norm(los);
                vs = obs.sat_velocities(:,sat_idx);
                doppler = obs.doppler(sat_idx);

                % Predicted velocity projection (unit: m/s)
                % pred = los' * R_total * v_s;
                % resi = -doppler + pred + bias;

                doppler_predicted =los'* (vs - R_total*local_vel) + bias;
                resi = doppler - doppler_predicted;

                % Construct Jacobian
                dRz_dyaw = [-sin(yaw), -cos(yaw), 0;
                             cos(yaw), -sin(yaw), 0;
                             0,         0,        0];
                dR_dyaw = R_rp * dRz_dyaw;

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
    R = R_rp * [cos(yaw), -sin(yaw), 0;
                sin(yaw),  cos(yaw), 0;
                0,         0,        1];
end
