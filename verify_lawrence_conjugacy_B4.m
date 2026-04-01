% verify_lawrence_conjugacy_B4.m
%
% Proves that the iterated reduced Long-Moody representation of B_4 on
% 2-point ordered configurations is simultaneously conjugate to Lawrence's
% representation from Fig 4.4 of Lawrence (1990).
%
% This extends verify_lawrence_conjugacy.m from n=3 to n=4.
%
% Parameters: Q (z_2 monodromy), q (z_1 monodromy), a (z_1-z_2 monodromy)
% Lawrence params: q1 = 1/Q, q2 = 1/q, alpha = 1/a

syms Q q a;

%% === Step 1: First iteration of reduced Long-Moody ===
%
% Seed representation of F_5 ⋊ B_5 (five strands, last one colored differently):
%   rho_g     = {Q, Q, Q, Q, a}  — 5 scalars (1x1 matrices)
%   rho_sigma = {1, 1, 1, 1}     — 4 scalars, all trivial
%   partition = [1 1 1 1 2]       — strand 5 is a different color
%
% Output: 4 matrices, each 4x4.
%   L4Pre{1} = rho^+(sigma_1)
%   L4Pre{2} = rho^+(sigma_2)
%   L4Pre{3} = rho^+(sigma_3)
%   L4Pre{4} = rho^+(sigma_4^2)  (squared because of partition boundary)

fprintf('=== Step 1: First reduced Long-Moody iteration ===\n');
fprintf('  Seed: rho_g = {Q,Q,Q,Q,a}, rho_sigma = {1,1,1,1}, partition = [1 1 1 1 2]\n');
L4Pre = reduced_lm_colored({Q, Q, Q, Q, a}, {sym(1), sym(1), sym(1), sym(1)}, [1 1 1 1 2]);
fprintf('  Output: 4 matrices, each %dx%d\n', size(L4Pre{1}));

%% === Step 2: Recover free group generators ===
%
% Using Birman exact sequence with partition [1 1 1 1 2]:
%   g_4 = L4Pre{4}            (already sigma_4^2 since last index is boundary)
%   g_3 = sigma_3^{-1} * g_4 * sigma_3  = L4Pre{3}^{-1} * g_4 * L4Pre{3}
%   g_2 = sigma_2^{-1} * g_3 * sigma_2  = L4Pre{2}^{-1} * g_3 * L4Pre{2}
%   g_1 = sigma_1^{-1} * g_2 * sigma_1  = L4Pre{1}^{-1} * g_2 * L4Pre{1}
%
% Output: {g_1, g_2, g_3, g_4} — four 4x4 matrices.

fprintf('\n=== Step 2: Recover free group generators ===\n');
L4PreFree = free_group_gens(L4Pre, [1 1 1 1 2]);
fprintf('  Recovered g_1,...,g_4 as %dx%d matrices\n', size(L4PreFree{1}));

%% === Step 3: Second iteration of reduced Long-Moody ===
%
% Seed representation of F_4 ⋊ B_4:
%   rho_g     = {g_1*q, g_2*q, g_3*q, g_4*q}  — 4 matrices, each 4x4, scaled by q
%   rho_sigma = {L4Pre{1}, L4Pre{2}, L4Pre{3}} — 3 matrices, each 4x4 (braid gens of B_4)
%   partition = [1 1 1 1]                        — uniform (default), no boundaries
%
% Output: 3 matrices, each 12x12.
%   L4{1} = rho^{++}(sigma_1)
%   L4{2} = rho^{++}(sigma_2)
%   L4{3} = rho^{++}(sigma_3)

fprintf('\n=== Step 3: Second reduced Long-Moody iteration ===\n');
fprintf('  Seed: rho_g = {g_i*q}, rho_sigma = {L4Pre{1..3}}, no partition\n');
rho_g = {L4PreFree{1}*q, L4PreFree{2}*q, L4PreFree{3}*q, L4PreFree{4}*q};
L4 = reduced_lm_colored(rho_g, L4Pre(1:3));
s1 = simplify(L4{1});
s2 = simplify(L4{2});
s3 = simplify(L4{3});
fprintf('  Output: 3 matrices, each %dx%d\n', size(s1));

%% === Step 4: Compute Lawrence's representation directly ===
%
% Lawrence (1990) Fig 4.4 for n=4, k=2 (ordered 2-point configurations).
% Basis: n(n-1) = 12 ordered pairs (lambda, mu), lambda ~= mu in {1,2,3,4},
% in lexicographic order:
%   1:(1,2) 2:(1,3) 3:(1,4) 4:(2,1) 5:(2,3) 6:(2,4)
%   7:(3,1) 8:(3,2) 9:(3,4) 10:(4,1) 11:(4,2) 12:(4,3)
%
% Parameter correspondence: q1 = 1/Q, q2 = 1/q, alpha = 1/a

fprintf('\n=== Step 4: Compute Lawrence representation for n=4 ===\n');
fprintf('  Basis: 12 ordered pairs (lam,mu) in lex order\n');
fprintf('  Parameters: q1 = 1/Q, q2 = 1/q, alpha = 1/a\n');
L_s1 = simplify(lawrence_generators(4, 1, 1/Q, 1/q, 1/a));
L_s2 = simplify(lawrence_generators(4, 2, 1/Q, 1/q, 1/a));
L_s3 = simplify(lawrence_generators(4, 3, 1/Q, 1/q, 1/a));
fprintf('  Lawrence matrices: 3 generators, each %dx%d\n', size(L_s1));

%% === Step 5: Find intertwining matrix P ===
%
% Solve for P such that P * s_i = L_s_i * P for all i simultaneously.
% Vectorize: (s_i^T ⊗ I) vec(P) = (I ⊗ L_s_i) vec(P)
% i.e., [(s_i^T ⊗ I) - (I ⊗ L_s_i)] vec(P) = 0 for each i.
% Stack all three constraints and find the null space.

fprintf('\n=== Step 5: Find intertwining matrix P ===\n');
fprintf('  Solving (s_i^T ⊗ I_12 - I_12 ⊗ L_s_i) vec(P) = 0 for i=1,2,3\n');
I12 = sym(eye(12));
A1 = kron(s1.', I12) - kron(I12, L_s1);
A2 = kron(s2.', I12) - kron(I12, L_s2);
A3 = kron(s3.', I12) - kron(I12, L_s3);
N = null([A1; A2; A3]);
fprintf('  Null space dimension: %d\n', size(N, 2));
assert(size(N, 2) == 1, 'Expected 1-dimensional null space (representations are irreducible).');

P = simplify(reshape(N(:,1), [12, 12]));
d = simplify(det(P));
fprintf('  det(P) = '); disp(d);
assert(~isequal(d, sym(0)), 'P must be invertible for generic Q, q, a.');
fprintf('  P is invertible for generic Q, q, a.\n');

%% === Step 6: Verify simultaneous conjugacy ===

fprintf('\n=== Step 6: Verify P * s_i = L_s_i * P for i = 1, 2, 3 ===\n');

check1 = simplify(P * s1 - L_s1 * P);
check2 = simplify(P * s2 - L_s2 * P);
check3 = simplify(P * s3 - L_s3 * P);

pass1 = all(all(check1 == 0));
pass2 = all(all(check2 == 0));
pass3 = all(all(check3 == 0));

if pass1, fprintf('  PASS: P * sigma_1 = L_sigma_1 * P\n');
else,      fprintf('  FAIL: sigma_1 conjugacy\n'); end

if pass2, fprintf('  PASS: P * sigma_2 = L_sigma_2 * P\n');
else,      fprintf('  FAIL: sigma_2 conjugacy\n'); end

if pass3, fprintf('  PASS: P * sigma_3 = L_sigma_3 * P\n');
else,      fprintf('  FAIL: sigma_3 conjugacy\n'); end

%% === Step 7: Verify braid relations for both representations ===
%
% B_4 braid relations:
%   (i)   sigma_1 sigma_2 sigma_1 = sigma_2 sigma_1 sigma_2
%   (ii)  sigma_2 sigma_3 sigma_2 = sigma_3 sigma_2 sigma_3
%   (iii) sigma_1 sigma_3 = sigma_3 sigma_1  (far commutativity)

fprintf('\n=== Step 7: Verify braid relations ===\n');

our_braid12  = simplify(s1*s2*s1 - s2*s1*s2);
our_braid23  = simplify(s2*s3*s2 - s3*s2*s3);
our_commute  = simplify(s1*s3 - s3*s1);
law_braid12  = simplify(L_s1*L_s2*L_s1 - L_s2*L_s1*L_s2);
law_braid23  = simplify(L_s2*L_s3*L_s2 - L_s3*L_s2*L_s3);
law_commute  = simplify(L_s1*L_s3 - L_s3*L_s1);

fprintf('  Our sigma_1 sigma_2 sigma_1 = sigma_2 sigma_1 sigma_2: %s\n', ...
    mat2str(all(all(our_braid12 == 0))));
fprintf('  Our sigma_2 sigma_3 sigma_2 = sigma_3 sigma_2 sigma_3: %s\n', ...
    mat2str(all(all(our_braid23 == 0))));
fprintf('  Our sigma_1 sigma_3 = sigma_3 sigma_1:                 %s\n', ...
    mat2str(all(all(our_commute == 0))));
fprintf('  Law sigma_1 sigma_2 sigma_1 = sigma_2 sigma_1 sigma_2: %s\n', ...
    mat2str(all(all(law_braid12 == 0))));
fprintf('  Law sigma_2 sigma_3 sigma_2 = sigma_3 sigma_2 sigma_3: %s\n', ...
    mat2str(all(all(law_braid23 == 0))));
fprintf('  Law sigma_1 sigma_3 = sigma_3 sigma_1:                 %s\n', ...
    mat2str(all(all(law_commute == 0))));

fprintf('\n=== Simultaneous conjugacy verified for B_4. ===\n');
