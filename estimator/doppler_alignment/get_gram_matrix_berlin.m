function M = get_gram_matrix_berlin(data)

n = size(data.gps_week, 1);
d = 3; % (TODO) only support 3 dimension now --sbs

M_obs = [];  
M = zeros(d^2 + 1 + 1, d^2 + 1 + 1);
for idx=1:n

    if mod(idx,10)~=0
       continue;
    end

    ecef_pos = data.ecef_position(idx, :)';
    local_vel = data.local_velocity(:,idx);

    sat_thres =4;
      max_sat = 0;
      if(size(data.sat_position{idx}, 2)<sat_thres)
          max_sat = size(data.sat_position{idx}, 2);
      else
          max_sat = sat_thres;
      end

     for sat = 1:max_sat
        % sat_idx = data{idx}.visible_sats(sat);
        sat_position = data.sat_position{idx}(:, sat);
        sat_velocity = data.sat_velocity{idx}(:, sat);
        doppler = -data.doppler{idx}(sat);
            if isnan(doppler)
                continue;
            end

        los_unit = (ecef_pos - sat_position)./norm(ecef_pos - sat_position,2);
        D = doppler + sat_velocity.'*los_unit;
        k = kron(local_vel', eye(d))'*los_unit;  
        M_obs = [M_obs; k'];  
        Mtt = get_mtt();
        Mrr = get_mrr(D, k);
        Mrt = get_mrt(D, k);
        M = M + [Mtt Mrt.'; Mrt Mrr];
    end

end


[U,S,~] = svd(M_obs);  
sing_vals = diag(S);
rank_obs = rank(M_obs); 
cond_num = max(sing_vals)/min(sing_vals);

rank_info.rank = rank_obs;
rank_info.singular_values = sing_vals;
rank_info.condition_number = cond_num;

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