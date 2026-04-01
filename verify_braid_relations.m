function ok = verify_braid_relations(rho_sigma, rho_g)
% VERIFY_BRAID_RELATIONS  Check braid relations and semidirect product relations.
%
%   ok = verify_braid_relations(rho_sigma)
%   ok = verify_braid_relations(rho_sigma, rho_g)
%
%   Given a cell array of matrices {M_1, ..., M_{n-1}} (candidate images of
%   the braid generators sigma_1, ..., sigma_{n-1}), this function checks:
%
%     1. Far commutativity:  M_i * M_j = M_j * M_i   for |i - j| > 1
%     2. Braid relation:    M_i * M_{i+1} * M_i = M_{i+1} * M_i * M_{i+1}
%
%   If rho_g is also provided (cell array of n matrices for the free group
%   generators g_1, ..., g_n), additionally checks the Bigelow-Tian
%   semidirect product relations of F_n ⋊ B_n:
%
%     3. g_{i+1} * sigma_i = sigma_i * g_i
%     4. g_i * sigma_i = sigma_i * g_i * g_{i+1} * g_i^{-1}
%     5. g_j * sigma_i = sigma_i * g_j   for j not in {i, i+1}
%
%   Returns true if all relations hold. Works for both symbolic and numeric
%   matrices.

    if ~iscell(rho_sigma)
        rho_sigma = num2cell(rho_sigma);
    end

    num_sigmas = length(rho_sigma);
    ok = true;

    is_symbolic = isa(rho_sigma{1}, 'sym');

    % --- Check far commutativity: sigma_i sigma_j = sigma_j sigma_i, |i-j|>1 ---
    for i = 1:num_sigmas
        for j = i+2:num_sigmas
            diff = rho_sigma{i} * rho_sigma{j} - rho_sigma{j} * rho_sigma{i};
            if ~check_zero(diff, is_symbolic)
                fprintf('FAIL: sigma_%d and sigma_%d do not commute.\n', i, j);
                ok = false;
            end
        end
    end

    % --- Check braid relation: sigma_i sigma_{i+1} sigma_i = sigma_{i+1} sigma_i sigma_{i+1} ---
    for i = 1:(num_sigmas - 1)
        lhs = rho_sigma{i} * rho_sigma{i+1} * rho_sigma{i};
        rhs = rho_sigma{i+1} * rho_sigma{i} * rho_sigma{i+1};
        if ~check_zero(lhs - rhs, is_symbolic)
            fprintf('FAIL: braid relation for sigma_%d, sigma_%d.\n', i, i+1);
            ok = false;
        end
    end

    if ok
        fprintf('All braid relations verified.\n');
    end

    % --- If free group generators provided, check semidirect product relations ---
    if nargin >= 2
        if ~iscell(rho_g)
            rho_g = num2cell(rho_g);
        end

        n = length(rho_g);
        % The semidirect product F_n ⋊ B_n uses sigma_1, ..., sigma_{n-1},
        % which are the first n-1 of the B_{n+1} generators.
        assert(num_sigmas >= n - 1, ...
            'Need at least %d braid generators for %d free group generators.', n-1, n);

        sdp_ok = true;

        for i = 1:(n-1)
            si = rho_sigma{i};
            gi = rho_g{i};
            gi1 = rho_g{i+1};
            gi_inv = inv(gi);

            % Relation 1: g_{i+1} * sigma_i = sigma_i * g_i
            if ~check_zero(gi1 * si - si * gi, is_symbolic)
                fprintf('FAIL: g_%d * sigma_%d ≠ sigma_%d * g_%d\n', i+1, i, i, i);
                sdp_ok = false;
            end

            % Relation 2: g_i * sigma_i = sigma_i * g_i * g_{i+1} * g_i^{-1}
            lhs = gi * si;
            rhs = si * gi * gi1 * gi_inv;
            if ~check_zero(lhs - rhs, is_symbolic)
                fprintf('FAIL: g_%d * sigma_%d ≠ sigma_%d * g_%d * g_%d * g_%d^{-1}\n', ...
                    i, i, i, i, i+1, i);
                sdp_ok = false;
            end

            % Relation 3: g_j * sigma_i = sigma_i * g_j for j not in {i, i+1}
            for j = 1:n
                if j == i || j == i+1
                    continue;
                end
                if ~check_zero(rho_g{j} * si - si * rho_g{j}, is_symbolic)
                    fprintf('FAIL: g_%d * sigma_%d ≠ sigma_%d * g_%d\n', j, i, i, j);
                    sdp_ok = false;
                end
            end
        end

        if sdp_ok
            fprintf('All semidirect product relations verified.\n');
        end
        ok = ok && sdp_ok;
    end
end


function z = check_zero(M, is_symbolic)
% CHECK_ZERO  Returns true if every entry of M is zero.
    if is_symbolic
        z = all(all(simplify(M) == 0));
    else
        z = (norm(double(M), 'fro') < 1e-10);
    end
end
