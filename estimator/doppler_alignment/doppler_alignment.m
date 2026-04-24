function [R,t, status, rank_info] = doppler_alignment(data, orthog, handed)
%ec_calibration Egomotion calibration.
status = 0;
if nargin < 2
    orthog = false;
end
if nargin < 3
    handed = false;
end

[M, rank_info] = get_gram_matrix(data);
% Get Schur complement
Mrr = M(2:end, 2:end);
Mtt = M(1, 1);
Mrt = M(2:end, 1);
Q = Mrr - Mrt*inv(Mtt)*Mrt.';

eigvals = eig(M);
min_eig = min(eigvals);


eigvals = eig(Q);
min_eig = min(eigvals);

[~, dual_sol] = dual_solver(Q, orthog, handed);

[R,status, rank_info.SVR] = extract_primal(dual_sol);
t = ec_get_t_opt(R, Mtt, Mrt.');
end
