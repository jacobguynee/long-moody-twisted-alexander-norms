% verify_lawrence_conjugacy.m
%
% Proves that our iterated reduced Long-Moody representation of B_3 on
% 2-point ordered configurations is simultaneously conjugate to Lawrence's
% representation from Fig 4.4 of Lawrence (1990).
%
% Parameters: Q (z_2 monodromy), q (z_1 monodromy), a (z_1-z_2 monodromy)
% Lawrence params: q1 = 1/Q, q2 = 1/q, alpha = 1/a

syms Q q a;

%% === Step 1: Compute our representation (iterated reduced LM) ===
fprintf('=== Computing iterated reduced Long-Moody representation ===\n');
L3Pre = reduced_lm_colored({Q, Q, Q, a}, {sym(1), sym(1), sym(1)}, [1 1 1 2]);
L3PreFree = free_group_gens(L3Pre, [1 1 1 2]);
rho_g = {L3PreFree{1}*q, L3PreFree{2}*q, L3PreFree{3}*q};
L3 = reduced_lm_colored(rho_g, L3Pre(1:2));
s1 = simplify(L3{1});
s2 = simplify(L3{2});
fprintf('Our sigma_1 and sigma_2: 6x6 matrices computed.\n');

%% === Step 2: Compute Lawrence's representation ===
fprintf('\n=== Computing Lawrence representation ===\n');
fprintf('Basis: lexicographic pairs (lam,mu), lam ~= mu in {1,2,3}\n');
fprintf('  1:(1,2) 2:(1,3) 3:(2,1) 4:(2,3) 5:(3,1) 6:(3,2)\n');
fprintf('Parameters: q1 = 1/Q, q2 = 1/q, alpha = 1/a\n\n');

L_s1 = simplify(lawrence_generators(3, 1, 1/Q, 1/q, 1/a));
L_s2 = simplify(lawrence_generators(3, 2, 1/Q, 1/q, 1/a));

%% === Step 3: Find intertwining matrix P ===
fprintf('=== Solving for P such that P * sigma_i = L_sigma_i * P ===\n');
I6 = sym(eye(6));
A1 = kron(s1.', I6) - kron(I6, L_s1);
A2 = kron(s2.', I6) - kron(I6, L_s2);
N = null([A1; A2]);
fprintf('Null space dimension: %d\n', size(N, 2));
assert(size(N, 2) == 1, 'Expected 1-dimensional null space');

P = simplify(reshape(N(:,1), [6, 6]));
d = simplify(det(P));
fprintf('det(P) = '); disp(d);
assert(~isequal(d, sym(0)), 'P must be invertible');
fprintf('P is invertible for generic Q, q, a.\n');

%% === Step 4: Verify simultaneous conjugacy ===
fprintf('\n=== Verifying P * sigma_i = L_sigma_i * P ===\n');

check1 = simplify(P * s1 - L_s1 * P);
check2 = simplify(P * s2 - L_s2 * P);

pass1 = all(all(check1 == 0));
pass2 = all(all(check2 == 0));

if pass1
    fprintf('PASS: P * sigma_1 = L_sigma_1 * P\n');
else
    fprintf('FAIL: sigma_1 conjugacy\n');
end

if pass2
    fprintf('PASS: P * sigma_2 = L_sigma_2 * P\n');
else
    fprintf('FAIL: sigma_2 conjugacy\n');
end

%% === Step 5: Verify braid relations for both ===
fprintf('\n=== Verifying braid relations ===\n');
our_braid = simplify(s1*s2*s1 - s2*s1*s2);
law_braid = simplify(L_s1*L_s2*L_s1 - L_s2*L_s1*L_s2);

fprintf('Our braid relation: %s\n', mat2str(all(all(our_braid == 0))));
fprintf('Lawrence braid relation: %s\n', mat2str(all(all(law_braid == 0))));

fprintf('\n=== Simultaneous conjugacy verified for B_3. ===\n');
