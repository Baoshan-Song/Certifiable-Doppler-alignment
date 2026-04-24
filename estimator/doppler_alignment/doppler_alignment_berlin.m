function [R,t, status] = doppler_alignment_berlin(data)
%ec_calibration Egomotion calibration.
status = 0;

M = get_gram_matrix_berlin(data);
% Get Schur complement
Mrr = M(2:end, 2:end);
Mtt = M(1, 1);
Mrt = M(2:end, 1);
Q = Mrr - Mrt*inv(Mtt)*Mrt.';

eigvals = eig(M);
min_eig = min(eigvals);


eigvals = eig(Q);
min_eig = min(eigvals);


% Use maximal linearly indpendent SO(3) constraints
orthog = true;
handed = true;

[~, dual_sol] = dual_solver(Q, orthog, handed);
[R,status] = extract_primal(dual_sol);
t = ec_get_t_opt(R, Mtt, Mrt.');
end
