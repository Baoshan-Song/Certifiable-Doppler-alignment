function H = da_hessian(muu1, gam, muu2, lam, muu3, lam_only)
% ec_hessian with selection constraints for rotation about z-axis
% muu1, muu2: multipliers for orthog constraints
% lam: multiplier for handedness
% muu3: multipliers for selection constraints (length 5)
% lam_only: flag to only compute lam part

d = size(muu1,1);

if nargin < 6
    lam_only = false;
end

% Original orthogonality Hessians
H_orthog1 = orthog_constraint_hessian(muu1);
H = [H_orthog1 zeros(d^2, 1); zeros(1, d^2) -trace(muu1)-gam];

%% Include X.'*X (duality strengthening)
if nargin > 2 && ~lam_only
    H_orthog2 = orthog_constraint_hessian(muu2, false);
    H = H + [H_orthog2 zeros(d^2, 1); zeros(1, d^2) -trace(muu2)];
end

%% Handedness constraint
if nargin > 3
    H_handed = handed_constraint_hessian(lam);
    H = H + H_handed;
end

%% Selection constraints (around z-axis)
if nargin > 4 && ~isempty(muu3)
    % vec(R) indexing (column-major, Matlab)
    % R13, R23, R31, R32, R33
    B = zeros(d^2 + 1, d^2 + 1, 5);
    idx = @(i,j) (j-1)*d + i; % vec(R) index

    n = d^2 + 1; % last row/col for lifting

    % 1. R13 = 0
    B(idx(1,3), n,1) = 1; B(n, idx(1,3),1) = 1;
    % 2. R23 = 0
    B(idx(2,3), n,2) = 1; B(n, idx(2,3),2) = 1;
    % 3. R31 = 0
    B(idx(3,1), n,3) = 1; B(n, idx(3,1),3) = 1;
    % 4. R32 = 0
    B(idx(3,2), n,4) = 1; B(n, idx(3,2),4) = 1;
    % 5. R33 = 1 (diagonal, only once)
    B(idx(3,3), n,5) = 1; % only fill one side, inner product gives M(idx,n)+M(n,idx)=2*M(idx,n)
    B(n, idx(3,3),5) = 1; % if you want symmetric

    % Add selection constraints with multipliers
    for q = 1:5
        H = H + muu3(q)*B(:,:,q);
    end
end

end
