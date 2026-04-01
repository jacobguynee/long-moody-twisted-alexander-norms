function rho_plus = reduced_lm_colored(rho_g, rho_sigma, partition)
% REDUCED_LM_COLORED  Reduced Long-Moody construction for colored/mixed braids.
%
%   rho_plus = reduced_lm_colored(rho_g, rho_sigma)
%   rho_plus = reduced_lm_colored(rho_g, rho_sigma, partition)
%
%   Like reduced_long_moody, but handles the case where the free group
%   generators g_i have different abelian labels (a "colored" coefficient
%   system). The representation rho is only well-defined on the mixed braid
%   group determined by the partition.
%
%   At a partition boundary (where strands i and i+1 are in different parts),
%   sigma_i is not in the valid subgroup, but sigma_i^2 is. The output for
%   such a boundary index is rho^+(sigma_i^2), computed correctly at the
%   group ring level:
%
%       phi_r(sigma_i^2) = sigma_i^2 * N_i' * N_i
%
%   where N_i' is N_i with g_i replaced by sigma_i^{-1}(g_i) = g_i g_{i+1} g_i^{-1}.
%
%   INPUTS:
%     rho_g     — cell array of n d×d matrices: {rho(g_1), ..., rho(g_n)}
%     rho_sigma — cell array of (n-1) d×d matrices: {rho(sigma_1), ..., rho(sigma_{n-1})}
%     partition — (optional) length-n vector of positive integers assigning
%                 each strand to a part, e.g. [1 1 1 2].
%                 Default: all ones (uniform partition, standard Long-Moody).
%
%   OUTPUT:
%     rho_plus  — cell array of (n-1) (n-1)*d × (n-1)*d matrices.
%                 rho_plus{i} = rho^+(sigma_i)   if i is not a partition boundary
%                 rho_plus{i} = rho^+(sigma_i^2) if i is a partition boundary
%
%   A partition boundary is an index i where partition(i) ~= partition(i+1).

    % --- Convert plain arrays to cell arrays for convenience ---
    if ~iscell(rho_g)
        rho_g = num2cell(rho_g);
    end
    if ~iscell(rho_sigma)
        rho_sigma = num2cell(rho_sigma);
    end

    n = length(rho_g);
    num_sigmas = length(rho_sigma);
    d = size(rho_g{1}, 1);
    N = (n - 1) * d;

    % --- Default partition: all strands in one group ---
    if nargin < 3 || isempty(partition)
        partition = ones(1, n);
    end
    assert(length(partition) == n, ...
        'Partition length (%d) must equal number of free group generators (%d).', ...
        length(partition), n);

    % --- Identify partition boundaries ---
    is_boundary = false(1, num_sigmas);
    for i = 1:num_sigmas
        is_boundary(i) = (partition(i) ~= partition(i+1));
    end

    % --- Validate inputs ---
    assert(num_sigmas == n - 1, ...
        'Expected %d braid generators (n-1), got %d.', n-1, num_sigmas);

    % --- Determine if we are working symbolically or numerically ---
    is_symbolic = isa(rho_g{1}, 'sym');
    if is_symbolic
        Id = sym(eye(d));
        Zd = sym(zeros(d));
    else
        Id = eye(d);
        Zd = zeros(d);
    end

    rho_plus = cell(1, num_sigmas);

    for i = 1:num_sigmas
        if ~is_boundary(i)
            % --- Interior of a partition part: standard sigma_i ---
            blocks = build_Ni(i, rho_g{i}, n, d, Id, Zd);

            si = rho_sigma{i};
            blocks = left_multiply_blocks(blocks, si, n-1);

            rho_plus{i} = assemble(blocks, n-1, d, N, is_symbolic);
        else
            % --- Partition boundary: compute sigma_i^2 at group ring level ---
            %
            % phi_r(sigma_i^2) = sigma_i^2 * N_i' * N_i
            %
            % N_i  uses g_i          in the S_i block.
            % N_i' uses g_i g_{i+1} g_i^{-1} in the S_i block.
            %
            % For d=1 (scalars commute): g_i g_{i+1} g_i^{-1} = g_{i+1}.
            % For d>1: rho(g_i g_{i+1} g_i^{-1}) = rho(g_i) * rho(g_{i+1}) * rho(g_i)^{-1}.

            gi_conj = rho_g{i} * rho_g{i+1} * inv(rho_g{i});  % rho(g_i g_{i+1} g_i^{-1})

            blocks_N  = build_Ni(i, rho_g{i}, n, d, Id, Zd);
            blocks_Np = build_Ni(i, gi_conj,  n, d, Id, Zd);

            % Multiply N_i' * N_i as block matrices.
            blocks_prod = multiply_block_matrices(blocks_Np, blocks_N, n-1, is_symbolic, d);

            % Left-multiply every block by rho(sigma_i)^2.
            si_sq = rho_sigma{i} * rho_sigma{i};
            blocks_prod = left_multiply_blocks(blocks_prod, si_sq, n-1);

            rho_plus{i} = assemble(blocks_prod, n-1, d, N, is_symbolic);
        end
    end
end


%% ========== Helper functions ==========

function blocks = build_Ni(i, gi, n, d, Id, Zd)
% BUILD_NI  Build the (n-1)×(n-1) block matrix N_i for the reduced S_i block.
%   gi is the d×d matrix to use in the S_i entries (either rho(g_i) or
%   the conjugated rho(g_i g_{i+1} g_i^{-1})).

    blocks = cell(n-1, n-1);
    for r = 1:(n-1)
        for c = 1:(n-1)
            if r == c
                blocks{r,c} = Id;
            else
                blocks{r,c} = Zd;
            end
        end
    end

    if i == 1
        % Boundary: sigma_1.
        %   [ -g   0 ]
        %   [  I   I ]
        blocks{1,1} = -gi;
        blocks{1,2} = Zd;
        blocks{2,1} = Id;
        blocks{2,2} = Id;

    elseif i == n - 1
        % Boundary: sigma_{n-1}.
        %   [  I   g ]
        %   [  0  -g ]
        blocks{n-2, n-2} = Id;
        blocks{n-2, n-1} = gi;
        blocks{n-1, n-2} = Zd;
        blocks{n-1, n-1} = -gi;

    else
        % Interior: sigma_i, 2 <= i <= n-2.
        %   [ I    g    0 ]
        %   [ 0   -g    0 ]
        %   [ 0    I    I ]
        p = i - 1;

        blocks{p,   p}   = Id;
        blocks{p,   p+1} = gi;
        blocks{p,   p+2} = Zd;

        blocks{p+1, p}   = Zd;
        blocks{p+1, p+1} = -gi;
        blocks{p+1, p+2} = Zd;

        blocks{p+2, p}   = Zd;
        blocks{p+2, p+1} = Id;
        blocks{p+2, p+2} = Id;
    end
end


function blocks = left_multiply_blocks(blocks, M, num_blocks)
% LEFT_MULTIPLY_BLOCKS  Left-multiply every block by M.
    for r = 1:num_blocks
        for c = 1:num_blocks
            blocks{r,c} = M * blocks{r,c};
        end
    end
end


function C = multiply_block_matrices(A, B, num_blocks, is_symbolic, d)
% MULTIPLY_BLOCK_MATRICES  Multiply two block matrices stored as cell arrays.
%   C{r,c} = sum_k A{r,k} * B{k,c}

    C = cell(num_blocks, num_blocks);
    for r = 1:num_blocks
        for c = 1:num_blocks
            if is_symbolic
                acc = sym(zeros(d));
            else
                acc = zeros(d);
            end
            for k = 1:num_blocks
                acc = acc + A{r,k} * B{k,c};
            end
            C{r,c} = acc;
        end
    end
end


function M = assemble(blocks, num_blocks, d, N, is_symbolic)
% ASSEMBLE  Assemble a cell array of blocks into one N×N matrix.
    if is_symbolic
        M = sym(zeros(N));
    else
        M = zeros(N);
    end
    for r = 1:num_blocks
        for c = 1:num_blocks
            rows = (r-1)*d+1 : r*d;
            cols = (c-1)*d+1 : c*d;
            M(rows, cols) = blocks{r,c};
        end
    end
    M = simplify(M);
end
