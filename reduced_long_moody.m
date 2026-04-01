function rho_plus = reduced_long_moody(rho_g, rho_sigma)
% REDUCED_LONG_MOODY  Reduced Long-Moody construction (Bigelow-Tian 2008, Section 6).
%
%   rho_plus = reduced_long_moody(rho_g, rho_sigma)
%
%   Given a representation rho of the semidirect product F_n ⋊ B_n,
%   encoded as matrices for the generators of F_n and B_n, this function
%   produces the reduced Long-Moody representation rho^+ of B_n.
%
%   INPUTS:
%     rho_g     — cell array of n d×d matrices: {rho(g_1), ..., rho(g_n)}
%                 (images of free group generators under rho)
%     rho_sigma — cell array of (n-1) d×d matrices: {rho(sigma_1), ..., rho(sigma_{n-1})}
%                 (images of braid group generators under rho)
%
%   OUTPUT:
%     rho_plus  — cell array of (n-1) (n-1)*d × (n-1)*d matrices:
%                 {rho^+(sigma_1), ..., rho^+(sigma_{n-1})}
%
%   The construction uses the 3×3 matrix S_i (with d×d block entries):
%
%       S_i = [ I    g_i    0 ]
%             [ 0   -g_i    0 ]
%             [ 0    I      I ]
%
%   placed at block positions (i-2, i-1, i) in an (n-1)×(n-1) block grid,
%   with identity blocks on the diagonal elsewhere. Boundary generators
%   sigma_1 and sigma_{n-1} use truncated versions of S_i. The entire
%   block matrix is then left-multiplied blockwise by rho(sigma_i).

    % --- Convert plain arrays to cell arrays for convenience ---
    %   Allows calling e.g. reduced_long_moody([q,q,q], [1,1]) for d=1.
    if ~iscell(rho_g)
        rho_g = num2cell(rho_g);
    end
    if ~iscell(rho_sigma)
        rho_sigma = num2cell(rho_sigma);
    end

    n = length(rho_g);           % number of free group generators
    num_sigmas = length(rho_sigma);  % should be n-1
    d = size(rho_g{1}, 1);      % dimension of the input representation
    N = (n - 1) * d;            % dimension of the output representation

    % --- Validate inputs ---
    assert(num_sigmas == n - 1, ...
        'Expected %d braid generators (n-1), got %d.', n-1, num_sigmas);
    for k = 1:n
        assert(all(size(rho_g{k}) == [d, d]), ...
            'rho_g{%d} should be %d×%d.', k, d, d);
    end
    for k = 1:num_sigmas
        assert(all(size(rho_sigma{k}) == [d, d]), ...
            'rho_sigma{%d} should be %d×%d.', k, d, d);
    end

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
        % Build the (n-1)×(n-1) block matrix for sigma_i.
        % Block indices are 1-based: block rows/cols 1, ..., n-1.
        %
        % The S_i block occupies block positions (i-1, i, i+1) when they
        % exist. (This is the 1-based version of the 0-based (i-2, i-1, i)
        % from the paper.)

        % Start with identity on the block diagonal.
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

        % Overwrite the S_i block region.
        % In 1-based block coords, S_i spans rows/cols from max(1, i-1) to min(n-1, i+1).
        gi = rho_g{i};  % rho(g_i)

        if i == 1
            % Boundary: sigma_1.
            % Truncated S_1 is 2×2 at block positions (1, 2):
            %   [ -g_1   0 ]
            %   [  I     I ]
            blocks{1,1} = -gi;
            blocks{1,2} = Zd;
            blocks{2,1} = Id;
            blocks{2,2} = Id;

        elseif i == n - 1
            % Boundary: sigma_{n-1}.
            % Truncated S_{n-1} is 2×2 at block positions (n-2, n-1):
            %   [  I    g_{n-1} ]
            %   [  0   -g_{n-1} ]
            blocks{n-2, n-2} = Id;
            blocks{n-2, n-1} = gi;
            blocks{n-1, n-2} = Zd;
            blocks{n-1, n-1} = -gi;

        else
            % Interior: sigma_i, 2 <= i <= n-2.
            % Full 3×3 S_i at block positions (i-1, i, i+1):
            %   [ I    g_i    0 ]
            %   [ 0   -g_i    0 ]
            %   [ 0    I      I ]
            p = i - 1;  % first block row/col of S_i (1-based)

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

        % Left-multiply every block by rho(sigma_i).
        si = rho_sigma{i};
        for r = 1:(n-1)
            for c = 1:(n-1)
                blocks{r,c} = si * blocks{r,c};
            end
        end

        % Assemble blocks into one N×N matrix.
        if is_symbolic
            M = sym(zeros(N));
        else
            M = zeros(N);
        end
        for r = 1:(n-1)
            for c = 1:(n-1)
                rows = (r-1)*d+1 : r*d;
                cols = (c-1)*d+1 : c*d;
                M(rows, cols) = blocks{r,c};
            end
        end
        rho_plus{i} = simplify(M);
    end
end
