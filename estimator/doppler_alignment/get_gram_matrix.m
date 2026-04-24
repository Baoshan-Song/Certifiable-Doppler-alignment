function [M, rank_info] = get_gram_matrix(data)

n = length(data);
d = 3; % (TODO) only support 3 dimension now --sbs

M_obs = []; 
V_obs = [];
N_obs =[];
M = zeros(d^2 + 1 + 1, d^2 + 1 + 1);
for idx=1:n

    ecef_pos = data{idx}.ecef_position;
    local_vel = data{idx}.local_velocity;
    for sat = 1:size(data{idx}.visible_sats)
        sat_idx = data{idx}.visible_sats(sat);
        sat_position = data{idx}.sat_positions(:, sat_idx);
        sat_velocity = data{idx}.sat_velocities(:, sat_idx);
        doppler = data{idx}.doppler(sat_idx);

        los_unit = (ecef_pos - sat_position)./norm(ecef_pos - sat_position,2);
        D = doppler + sat_velocity.'*los_unit;
        k = kron(local_vel', eye(d))'*los_unit;
        M_obs = [M_obs; k']; 
        N_obs = [N_obs;los_unit'];
        V_obs = [V_obs;local_vel'];
        Mtt = get_mtt();
        Mrr = get_mrr(D, k);
        Mrt = get_mrt(D, k);
        M = M + [Mtt Mrt.'; Mrt Mrr];
    end

end

[U,S,~] = svd(M);  
tol = 2e-7;
sing_vals = diag(S);
real_rank = sum(sing_vals > tol);
rank_obs = rank(M);  
cond_num = max(sing_vals)/min(sing_vals);

rank_info.rank_M = real_rank;
rank_info.singular_values_M = sing_vals;
rank_info.condition_number_M = cond_num;


% [U,S,~] = svd(N_obs);  
% tol = 1e-6;
% sing_vals = diag(S);
% real_rank = sum(sing_vals > tol);
% rank_obs = rank(N_obs); 
% cond_num = max(sing_vals)/min(sing_vals);
% rank_info.rank_N = real_rank;
% rank_info.singular_values_N = sing_vals;
% rank_info.condition_number_N = cond_num;

% [U,S,~] = svd(V_obs);  
% tol = 1e-6;
% sing_vals = diag(S);
% real_rank = sum(sing_vals > tol);
% rank_obs = rank(V_obs); 
% cond_num = max(sing_vals)/min(sing_vals);
% rank_info.rank_V = real_rank;
% rank_info.singular_values_V = sing_vals;
% rank_info.condition_number_V = cond_num;

end

function Mtt = get_mtt()
Mtt = 1;
end

function Mrr = get_mrr(D, k)
Mrr = k*k';
Mrc = -D*k;
Mrr = [Mrr Mrc;Mrc.' D*D];
end

function Mrt = get_mrt(D, k)
Mrt = [k; -D];

end