% BURAU_EXAMPLE  Compute the reduced Burau representation via the reduced
%                Long-Moody construction and verify braid relations.
%
%   The reduced Burau arises from the simplest input: a one-dimensional
%   representation rho of F_n ⋊ B_n given by
%       rho(g_i) = q   (symbolic variable)
%       rho(sigma_i) = 1
%
%   By Bigelow-Tian Proposition 2.1 and Section 6, the Long-Moody
%   construction applied to this input yields the unreduced Burau, and
%   the reduced version yields the reduced Burau.

%% === Parameters ===
n = 4;  % number of strands in B_n (change this to try other sizes)

%% === Set up the input representation ===
syms q

% rho(g_i) = q for all i = 1, ..., n  (1×1 symbolic matrices)
rho_g = cell(1, n);
for i = 1:n
    rho_g{i} = q;
end

% rho(sigma_i) = 1 for all i = 1, ..., n-1
rho_sigma = cell(1, n-1);
for i = 1:(n-1)
    rho_sigma{i} = sym(1);
end

%% === Run the reduced Long-Moody construction ===
fprintf('Reduced Long-Moody construction for B_%d with Burau input\n', n);
fprintf('Input: rho(g_i) = q, rho(sigma_i) = 1, d = 1\n');
fprintf('Output: %d matrices of size %d x %d\n\n', n-1, n-1, n-1);

rho_plus = reduced_long_moody(rho_g, rho_sigma);

%% === Display the output matrices ===
for i = 1:(n-1)
    fprintf('rho^+(sigma_%d) =\n', i);
    disp(rho_plus{i});
end

%% === Verify braid relations ===
fprintf('Verifying braid relations...\n');
verify_braid_relations(rho_plus);
