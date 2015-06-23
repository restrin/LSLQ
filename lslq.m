function [x, flag, relres, iter, resvec] = lslq( A, b, atol, btol, conlim, maxit )
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
%   [X,FLAG,RELRES,ITER] = LSLQ(A,B,...) also returns the iteration
%   number at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = SYMMLQ(A,B,...) also returns a vector of
%   of estimates of the SYMMLQ residual norms at each iteration, including
%   NORM(B-A*X0).


end