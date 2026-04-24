close all; clear; clc;
rng(41);
simulation;


%% Global: Construct and solve semidefinite relaxation
disp('Certifiable estimation:');
tic;
[R, dt, status, rank_info] = doppler_alignment(data_cells);
% [R, dt, status, rank_info] = doppler_alignment(data_cells, true, true);
time1 = toc
if status == 0
    disp('Angular Error to Ground Truth (deg):');
    angle_err_rad = norm(acos((trace(R * R_ecef2local) - 1)/2));
    disp(rad2deg(angle_err_rad));
else
   disp('Rank tightness fail:'); 
    disp('Inexact SDP Angular Error to Ground Truth (deg):');
    angle_err_rad = norm(acos((trace(R * R_ecef2local) - 1)/2));
    disp(rad2deg(angle_err_rad));
end

%% Local 1: Solve rotation with raw Doppler 
disp('Gauss-Newton estimation:');
tic
[R_gn, b_gn] = gauss_newton(data_cells, lambda , eye(3));
time2 = toc
disp('Estimated R_e_w:');
disp(R_gn');
disp('Estimated Clock Drift (m/s):');
disp(b_gn);
disp('Angular Error to Ground Truth (deg):');
angle_err_rad = acos((trace(R_gn * R_ecef2local) - 1)/2);
disp(rad2deg(angle_err_rad));

%% Local 2: Solve rotation with SPV and Wahba 
disp('SPV + Velocity registration:');
tic;
[R_icp, b_icp] = estimate_rotation_from_doppler(data_cells);
time3 = toc
disp('Estimated R_e_w:');
disp(R_icp');
disp('Estimated Clock Drift (m/s):');
disp(b_icp);
disp('Angular Error to Ground Truth (deg):');
angle_err_rad = acos((trace(R_icp * R_ecef2local) - 1)/2);
disp(rad2deg(angle_err_rad));




function [R_est, b_est] = estimate_rotation_from_doppler(data_cells)
    V_local = [];
    V_est = [];
    b_all = [];

    for k = 1:length(data_cells)
        data = data_cells{k};
        vis = data.visible_sats;
        A = [];
        y = [];

        for j = 1:length(vis)
            v_r_local = data.local_velocity;  % current local vel
            
            s = vis(j);
            los = data.ecef_position - data.sat_positions(:, s);  % rx - sat
            los_unit = los / norm(los);

            vs = data.sat_velocities(:, s);  % v_s^e
            z = data.doppler(s);             % [m/s]

            A = [A; los_unit'];  % 对 R*v_r
            y = [y; z + los_unit' * vs];
        end

        if size(A,1) < 3
            continue;  
        end

        A_full = [A, ones(size(A,1),1)];
        sol = A_full \ y;

        x_est = sol(1:3);  % R * v_r_local
        b_est_k = sol(4);

        V_local = [V_local, v_r_local];
        V_est = [V_est, x_est];
        b_all = [b_all; b_est_k];
    end

    %  V_est ≈ R * V_local
    M = V_est * V_local';  
    [U, ~, V] = svd(M);
    R_est = U * V';

    if det(R_est) < 0
        U(:,end) = -U(:,end);
        R_est = U * V';
    end

    b_est = mean(b_all);  
end
