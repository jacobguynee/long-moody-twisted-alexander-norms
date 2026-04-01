function rho_g = free_group_gens(rho_sigma, partition)
% FREE_GROUP_GENS  Compute images of free group generators from braid generators.
%
%   rho_g = free_group_gens(rho_sigma)
%   rho_g = free_group_gens(rho_sigma, partition)
%
%   Given the output of reduced_long_moody or reduced_lm_colored (n matrices
%   for the generators of B_{n+1} or the mixed braid group), compute the
%   images of the n free group generators g_1, ..., g_n of F_n using the
%   recursive formula:
%
%       g_n = sigma_n^2
%       g_i = sigma_i^{-1} * g_{i+1} * sigma_i   (for i = n-1, ..., 1)
%
%   This is equivalent to the Birman exact sequence formula
%       g_i = (sigma_n ... sigma_{i+1}) * sigma_i^2 * (sigma_{i+1}^{-1} ... sigma_n^{-1})
%   but avoids using sigma_n directly in intermediate steps, so it works
%   correctly with colored output where the last entry may be sigma_n^2.
%
%   PARTITION:
%     Optional length-(n+1) vector assigning strands to partition parts,
%     e.g. [1 1 1 2]. When partition(n) ~= partition(n+1), the last entry
%     of rho_sigma is already sigma_n^2 (from reduced_lm_colored) and is
%     used directly as g_n. Otherwise, sigma_n is squared to get g_n.
%
%     NOTE: This function currently only handles partitions where the sole
%     boundary is at the last index — i.e., partitions of the form
%     [1 1 ... 1] (uniform) or [1 1 ... 1 2] (last strand pure).
%     General partitions with multiple interior boundaries would require
%     the actual sigma_i matrices at those boundaries, which the colored
%     output does not provide.
%
%   INPUTS:
%     rho_sigma — cell array (or plain array for d=1) of n d×d matrices.
%     partition — (optional) length-(n+1) vector, default all ones.
%
%   OUTPUT:
%     rho_g — cell array of n d×d matrices: {rho(g_1), ..., rho(g_n)}

    % Convert plain array to cell array for convenience.
    if ~iscell(rho_sigma)
        rho_sigma = num2cell(rho_sigma);
    end

    n = length(rho_sigma);  % number of generators

    % Default partition: uniform (all strands in one group).
    if nargin < 2 || isempty(partition)
        partition = ones(1, n + 1);
    end

    % Determine whether the last entry is already sigma_n^2.
    last_is_boundary = (partition(n) ~= partition(n + 1));

    % Base case: g_n = sigma_n^2.
    if last_is_boundary
        % The last entry of rho_sigma is already sigma_n^2.
        gn = rho_sigma{n};
    else
        % The last entry is sigma_n; square it.
        gn = rho_sigma{n} * rho_sigma{n};
    end

    rho_g = cell(1, n);
    rho_g{n} = simplify(gn);

    % Recursive step: g_i = sigma_i^{-1} * g_{i+1} * sigma_i.
    % (Right conjugation — matches the Artin action convention.)
    for i = (n-1):-1:1
        rho_g{i} = simplify(inv(rho_sigma{i}) * rho_g{i+1} * rho_sigma{i});
    end
end
