function M = lawrence_generators(n, i, q1, q2, alpha)
% LAWRENCE_GENERATORS  Lawrence's matrix for sigma_i on ordered 2-point configurations.
%
%   M = lawrence_sigma(n, i, q1, q2, alpha)
%
%   Returns the n(n-1) x n(n-1) matrix for the action of the braid generator
%   sigma_i on the homology of the ordered configuration space of 2 points
%   in the n-punctured disk, from Lawrence (1990) Fig 4.4.
%
%   Basis: ordered pairs (w_lambda, w_mu) with lambda ~= mu in {1,...,n},
%   in lexicographic order.
%
%   Parameters q1, q2 are the monodromies of the two configuration points
%   around punctures; alpha is the monodromy of one point around the other.

    d = n*(n-1);
    M = sym(zeros(d));

    % Map (lambda, mu) -> index in lexicographic order.
    idx = containers.Map('KeyType','char','ValueType','int32');
    pairs = zeros(d, 2);
    count = 0;
    for lam = 1:n
        for mu = 1:n
            if lam ~= mu
                count = count + 1;
                pairs(count,:) = [lam, mu];
                idx(sprintf('%d,%d', lam, mu)) = count;
            end
        end
    end

    % Helper to get index of pair (a,b).
    p = @(a,b) idx(sprintf('%d,%d', a, b));

    ip1 = i + 1;  % i+1

    % --- Row w_{i,i+1} ---
    r = p(i, ip1);
    M(r, p(ip1, i)) = 1/(q1*q2);
    for j = 1:(i-1)
        M(r, p(ip1, j)) = (1/q2 - 1)/q1;
    end
    for k = (ip1+1):n
        M(r, p(ip1, k)) = (1/q2 - 1)/q1;
    end

    % --- Row w_{i+1,i} ---
    r = p(ip1, i);
    M(r, p(i, ip1)) = 1/(q1*q2*alpha);
    for j = 1:(i-1)
        M(r, p(j, ip1)) = (1/q1 - 1)/(q2*alpha);
    end
    for k = (ip1+1):n
        M(r, p(k, ip1)) = (1/q1 - 1)/q2;
    end

    % --- For each spectator j < i ---
    for j = 1:(i-1)
        % Row w_{ij}: 0*w_{ij} + (1/q1)*w_{i+1,j}
        r = p(i, j);
        M(r, p(i, j)) = 0;  % explicit zero (default)
        M(r, p(ip1, j)) = 1/q1;

        % Row w_{ji}: 0*w_{ji} + (1/q2)*w_{j,i+1}
        r = p(j, i);
        M(r, p(j, i)) = 0;
        M(r, p(j, ip1)) = 1/q2;

        % Row w_{i+1,j}: 1*w_{ij} + (1-1/q1)*w_{i+1,j}
        r = p(ip1, j);
        M(r, p(i, j)) = 1;
        M(r, p(ip1, j)) = 1 - 1/q1;

        % Row w_{j,i+1}: 1*w_{ji} + (1-1/q2)*w_{j,i+1}
        r = p(j, ip1);
        M(r, p(j, i)) = 1;
        M(r, p(j, ip1)) = 1 - 1/q2;
    end

    % --- For each spectator k > i+1 ---
    for k = (ip1+1):n
        % Row w_{ik}: 0*w_{ik} + (1/q1)*w_{i+1,k}
        r = p(i, k);
        M(r, p(i, k)) = 0;
        M(r, p(ip1, k)) = 1/q1;

        % Row w_{ki}: 0*w_{ki} + (1/q2)*w_{k,i+1}
        r = p(k, i);
        M(r, p(k, i)) = 0;
        M(r, p(k, ip1)) = 1/q2;

        % Row w_{i+1,k}: 1*w_{ik} + (1-1/q1)*w_{i+1,k}
        r = p(ip1, k);
        M(r, p(i, k)) = 1;
        M(r, p(ip1, k)) = 1 - 1/q1;

        % Row w_{k,i+1}: 1*w_{ki} + (1-1/q2)*w_{k,i+1}
        r = p(k, ip1);
        M(r, p(k, i)) = 1;
        M(r, p(k, ip1)) = 1 - 1/q2;
    end

    % --- Diagonal 1's for pairs with both indices not in {i, i+1} ---
    for s = 1:d
        lam = pairs(s, 1);
        mu  = pairs(s, 2);
        if lam ~= i && lam ~= ip1 && mu ~= i && mu ~= ip1
            M(s, s) = 1;
        end
    end

    M = simplify(M);
end
