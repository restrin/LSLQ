function [ x, iter, resvec ] = lslqnaive3( A, b, tol, maxiter )

resvec = zeros(maxiter,1);

m = size(A,1);
n = size(A,2);

% W = [w ww]
W = zeros(n, 2);

x = zeros(n,1);

nb = norm(b);
u = b/nb;
v = A'*u;
alpha = norm(v);
v = v/alpha;
W(:,2) = v;

k = 1;

alpha1 = alpha;
alphap = alpha;
betap = nb;

while ( k < maxiter )
    
    u = A*v - alpha*u;
    beta = norm(u);
    u = u/beta;

    v = A'*u - beta*v;
    alpha = norm(v);
    v = v/alpha;

    if (k >= 2)
        % Estimating norm of A'*r
        Atr = [alphap*betap (alphap^2 + beta^2); 0 alpha*beta]*y;
        nAtr = norm(Atr);
        
 %       norm(A'*A*x - A'*b) - nAtr
        
        resvec(k) = nAtr;
        
        if (nAtr < tol*alpha1*nb)
            break;
        end
    end
    
    if (k == 1)        
        rho = sqrt(alphap^2 + beta^2);
        c1 = alphap / rho;
        s1 = -beta / rho;
        
        %%% FOR RESIDUAL
            q = [c1 s1];
        %%%
        
        sigmahat = rho;
        
        % theta == p
        theta = alpha * beta / rho;        
        z = alphap * nb / rho;
        
        sigma = sqrt(sigmahat^2 + theta^2);
        c2 = sigmahat / sigma;
        s2 = -theta / sigma;
        
        zz = z/sigma;
 
        W(:,1) = W(:,2);
        W(:,2) = v;
        W = W*[c2, s2; -s2, c2];
        x = zz*W(:,1);
        
        y = [c2; -s2] * zz;
    else              
        
        % Q
%        thetam1 = theta; %-s1*alphap;
        rhohat = c1*alphap;
        
        rho = sqrt(rhohat^2 + beta^2);
        c1 = rhohat/rho;
        s1 = -beta/rho;     
        
        %%%%%%%%%
        % RESIDUAL NORM CALCULATION for iteration k-1
            % Doesn't work for iteration k = 1
            q = q(2)*[-c1; s1];
            ry = [rho*y(2); 0];
            r = norm(q*nb - ry);
            
            norm(A*x - b) - r
        %%%%%%%%%
        
        theta = alpha*beta/rho;
        z = -thetam1/rho * z;
        
        eta = -s2*rho;
        sigmahat = c2*rho;
        
        % Qtilde
        sigma = sqrt(sigmahat^2 + theta^2);
        c2 = sigmahat/sigma;
        s2 = -theta/sigma;
        
        zz = (z - eta*zz)/sigma;
        
        W(:,1) = W(:,2);
        W(:,2) = v;
        W = W*[c2 s2; -s2 c2];
        
        x = x + zz*W(:,1);
        
        % To get norm of residual      
        y = [sp -cp*c2;0 s2]*[zzp; zz];
    end
    
    alphap = alpha;
    betap = beta;
    cp = c2;
    sp = s2;
    zzp = zz;
    thetam1 = theta;
    
    k = k+1;
end

iter = k;

end

