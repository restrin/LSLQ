function [x, flag, relres, normAr, iter, resvec] = lslq( A, b, atol, btol, conlim, maxit )
% LSLQ Least-Squares LQ method
%   X = LSLQ(A,B) attempts to solve the system of linear equations A*X=B
%   for X. B is a column vector of length N. The system must be consistent.
%
%   X = LSLQ(AFUN,B) accepts a function handle AFUN instead of the matrix
%   A. AFUN(X) accepts a vector input X and returns the matrix-vector
%   product A*X. In all of the following syntaxes, you can replace A by
%   AFUN.
%
%   X = LSLQ(A,B,ATOL,BTOL) continues iterations until a certain
%   backward error estimate is smaller than some quantity depending on
%   ATOL and BTOL.  Let RES = B - A*X be the residual vector for the
%   current approximate solution X.  If A*X = B seems to be consistent,
%   LSLQ terminates when NORM(RES) <= ATOL*NORM(A)*NORM(X) + BTOL*NORM(B).
%   If the system is inconsistent, we will throw some error. INCLUDE
%   DOCUMENTATION ON THIS.
% 
%   X = LSLQ(A,B,ATOL,BTOL,CONLIM) terminates if an estimate
%   of cond(A) exceeds CONLIM. For compatible systems Ax = b,
%   conlim could be as large as 1.0e+12 (say). If CONLIM is [], the default
%   value is CONLIM = 1e+8. Maximum precision can be obtained by setting 
%   ATOL = BTOL = CONLIM = 0, but the number of iterations may then be 
%   excessive.
%
%   X = LSLQ(A,B,ATOL,BTOL,MAXIT) specifies the maximum number of
%   iterations. If MAXIT is [] then SYMMLQ uses the default, min(N,20).
%
%   [X,FLAG] = LSLQ(A,B,...) also returns a convergence FLAG:
%    0 LSLQ converged to the desired tolerance
%      NORM(RES) <= ATOL*NORM(A)*NORM(X) + BTOL*NORM(B) within MAXIT
%      iterations.
%    1 LSLQ iterated MAXIT times but did not converge.
%    MORE FLAGS TO COME
%
%   [X,FLAG,RELRES] = LSLQ(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then
%   RELRES <= ATOL*NORM(A)*NORM(X) + BTOL*NORM(B).
%
%   [X,FLAG,RELRES, NORMAR] = LSLQ(A,B,...) also returns A'*R, where
%   R = B - A*X. If FLAG is 0, then
%   A'*RELRES <= NORM(A)*ATOL*NORM(A)*NORM(X) + NORM(A)*BTOL*NORM(B).
%
%   [X,FLAG,RELRES,ITER] = LSLQ(A,B,...) also returns the iteration
%   number at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = LSLQ(A,B,...) also returns a vector of
%   of estimates of the LSLQ residual norms at each iteration, including
%   NORM(B-A*X0).

% Initialization
  if isa(A,'numeric')
    explicitA = true;
  elseif isa(A,'function_handle')
    explicitA = false;
  else
    error('SOL:lsmr:Atype','%s','A must be numeric or a function handle');
  end
  
  % Determine dimensions m and n, and
  % form the first vectors u and v.
  % These satisfy  beta*u = b,  alpha*v = A'u.

  u    = b;
  beta = norm(u);
  n2b = beta;
  if beta > 0
    u  = u/beta;
  end

  if explicitA
    v = A'*u;
    [m n] = size(A);
  else  
    v = A(u,2);
    m = size(b,1);
    n = size(v,1);
  end

  minDim = min([m n]);
  
  % Set default parameters.
  
  if nargin < 3 || isempty(atol)     , atol      = 1e-6;       end
  if nargin < 4 || isempty(btol)     , btol      = 1e-6;       end
  if nargin < 5 || isempty(conlim)   , conlim    = 1e+8;       end
  if nargin < 6 || isempty(maxit)    , maxit     = minDim;     end
  
  alpha = norm(v);
  if alpha > 0
    v = (1/alpha)*v;
  end
  
  anorm = alpha;

  W = zeros(n, 2);
  W(:,1) = v;

  alphap = alpha;
%  betap = beta;
  
  % First iteration
  % Golub-Kahan step
  if (explicitA)
    u = A * v  - alpha*u;
  else
    u = A(v,1) - alpha*u; 
  end
  beta = norm(u);
  
  if (beta > 0)
    u = u/beta;

    if (explicitA)
      v = A' * u - beta*v;
    else
      v = A(u,2) - beta*v;
    end
    alpha = norm(v);
    if (alpha > 0)
      v = v/alpha;
    end
  end
    
  % First iteration
  rho = sqrt(alphap^2 + beta^2);
  c1  = alphap / rho;
  s1  = -beta / rho;

  q = [c1 s1];

  sigmahat = rho;

  theta = alpha * beta / rho;        
  z     = alphap * n2b / rho;

  sigma = sqrt(sigmahat^2 + theta^2);
  c2    = sigmahat / sigma;
  s2    = -theta / sigma;

  zz    = z/sigma;

  W(:,2) = v;
  W      = W*[c2, s2; -s2, c2];
  x      = zz*W(:,1);
  
  y = [c2; -s2] * zz;
  
  resvec = zeros(maxit,1);
  
  for it = 2:maxit
    
    alphap = alpha;
    betap = beta;
    cp = c2;
    sp = s2;
    zzp = zz;
    thetap = theta;
      
    % Golub-Kahan step
    if (explicitA)
      u = A * v  - alpha*u;
    else
      u = A(v,1) - alpha*u; 
    end
    beta = norm(u);
  
    if (beta > 0)
      u = u/beta;

      if (explicitA)
        v = A' * u - beta*v;
      else
        v = A(u,2) - beta*v;
      end
      alpha = norm(v);
      if (alpha > 0)
        v = v/alpha;
      end
    end
    
    anorm = sqrt(anorm^2 + alpha^2 + beta^2);

    % Estimating norm of A'*r
    Atr = [alphap*betap (alphap^2 + beta^2); 0 alpha*beta]*y;
    n2Atr = norm(Atr);
    
    rhohat = c1*alphap;

    rho = sqrt(rhohat^2 + beta^2);
    c1 = rhohat/rho;
    s1 = -beta/rho;     

    % Residual Estimation
    q = q(2)*[-c1; s1];
    ry = [rho*y(2); 0]*((-1)^(sign(it-2)-1)); % First iteration is different for some reason
    n2r = norm(q*n2b - ry);
    resvec(it-1) = n2r;
    
    if (n2r < atol*anorm + btol*n2b)
      flag   = 0;
      iter   = it-1;
      relres = n2r;
      normAr = n2Atr;
      resvec = resvec(1:it-1);
      return;
    end

    theta = alpha*beta/rho;
    z = -thetap/rho * z;

    eta = -s2*rho;
    sigmahat = c2*rho;

    sigma = sqrt(sigmahat^2 + theta^2);
    c2 = sigmahat/sigma;
    s2 = -theta/sigma;

    zz = (z - eta*zz)/sigma;

    W(:,1) = W(:,2);
    W(:,2) = v;
    W = W*[c2 s2; -s2 c2];

    x = x + zz*W(:,1);
     
    y = [sp -cp*c2;0 s2]*[zzp; zz];
    
  end
  
% We exceeded the maximum number of iterations
flag   = 1;
iter   = it;
relres = n2r;
normAr = n2Atr;
  
end