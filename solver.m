function [ydual, solverTime] = solver(mode, A, b, K, Q, tau)
%-------------------------------------------------------------
% Provide different solvers for user to solve dual problem
% mode 1: built-in quadprog
% mode 2: pdco by Michael Saunder
% mode 3: ASP by Michael P. Friedlander and Michael Saunder
% 
% update
% Mar 1: - pdco: d1 from sparse to zeros, because print option doesn't
%                support sparse vector. 
%                e.g. d1 = sparse(5,1); d1(4) =1; max(d1)= (4,1) 1;
%        - pdco: interior method cannot warm-start, thus delete v0=zeros
%-------------------------------------------------------------

    [m,n] = size(A);
    switch mode
        case 1
            % quadprob with "interior-point-convex" algorithm
            % min_x 1/2 xHx+f'x 
            % subject to Ax<=b or Aeq x = beq or lb <= x <= up
            
            options = optimoptions(@quadprog,'Display','off'); 
            H = [eye(K+2), zeros(K+2,1);zeros(1,K+2) 0]; 
            f = [Q'*b;-tau];
            A = [A'*Q -ones(n,1);-A'*Q -ones(n,1)]; 
            b = zeros(2*n,1); 
            
            quadTimeVal = tic();
            [v,fquad,exitflag] = quadprog(H,f,A,b,[],[],[],[],[],options); 
            solverTime = toc(quadTimeVal);

            ydual = Q*v(1:K+2,1);
            fprintf('quadprog time is : %16e \n', solverTime);
            
        case 2
            % pdco (Primal-Dual Barrier Method)
            %    minimize    phi(x) + 1/2 norm(D1*x)^2 + 1/2 norm(r)^2
            %      x,r
            %    subject to  A*x + D2*r = b,   bl <= x <= bu,   r unconstrained
            
           
            options = pdcoSet;
            options.Method = 22;
            options.Print = 0; % no print

            
            %H = sparse(K+2+1+2*n, K+2+1+2*n); 
            %H(1:K+2, 1:K+2) = speye(K+2); %= Q'*Q is an identity matrix
            %objectivefunction = @(x) deal(0.5*(x'*H*x) + c'*x, H*x + c, H);
            c = [-Q'*b; tau; sparse(2*n,1)];
            A = [-A'*Q, -ones(n,1), speye(n), sparse(n,n); A'*Q, -ones(n,1),sparse(n,n), speye(n)];
            b = sparse(2*n, 1); 
            bl = [-Inf(K+2+1,1);sparse(2*n,1)]; 
            bu = Inf(K+2+1+2*n,1);
            
            % create warm start use previous iteration solution
            
            v0 = [0;zeros(K+2+2*n,1)]; % warm start
            %A*v0 == b
            
            d1 = zeros(K+2+1+2*n,1); 
            d1(1:K+2) = 1;

            pdcoTimeVal = tic();
            %include quadratic part into the obj function
            %[v,y,z,inform,PDitns,CGitns,time]=pdco(objectivefunction,Amatrix,bvec,bl,bu,1e-4,1e-4,options1,v0,sparse(2*n,1),sparse(K+1+2*n,1),1,1);
            %include quadratic part into the regularization D1
            [v,y,z,inform,PDitns,CGitns,time] = pdco(c,A,b,bl,bu,d1,1e-4,options,v0,sparse(2*n,1),sparse(K+2+1+2*n,1),1,1);
            %keyboard
            solverTime   = toc(pdcoTimeVal);
            

            %print info about dual subproblem
            switch inform 
                case 0
                    fprintf('\n Dual solution is found \n');
                case 1
                    fprintf('\n Too many iterations were required, exceed maximum iterations \n ');
                case 2
                    fprintf('\n Linesearch failed too often \n ');
                case 3
                    fprintf('\n The step lengths became too small \n');
                otherwise
                    fprintf('\n Cholesky said ADDA was not positive definite \n ');
            end
            ydual = Q*v(1:K+2); %v=[c,lambda, s1, s2], c: (K+2)x 1
            
            %fprintf('pdco time is : %16e \n', solverTime);
            %fprintf('obj value is %15e \n',1/2*(ydual'*ydual)+c(1:K+2)'*v(1:K+2)+tau*v(K+3));
            
        case 3
            % ASP (active-set pursuit)
            % ---------------------------------
            % min_y 1/2 lambda y'y -b'y 
            % subject to bl <= A'y <= bu
            % ---------------------------------
            
            [ydual, solverTime, inform] = ASPwrap(A,b,Q,tau,K);
            ydual = Q*ydual;
            fprintf('ASP time is : %16e \n', solverTime);
            
        
        case 4
            % Gurobi
            
            Q1 = sparse(2*n+K+2+1,2*n+K+2+1);
            Q1(1:K+2,1:K+2) = 1/2.*speye(K+2);

            model.Q = Q1;
            model.obj = [-Q'*b; tau; zeros(2*n,1)];
            model.A = sparse([-A'*Q -ones(n,1) eye(n) zeros(n); A'*Q -ones(n,1) zeros(n) eye(n)]);
            model.rhs = zeros(2*n,1);
            model.sense = '=';
            lb = zeros(2*n+K+2+1,1); 
            lb(1:K+2+1) = -inf; 
            model.lb = lb;
            gurobi_write(model, 'qp.lp');
            params.outputflag = 0; 
            results = gurobi(model,params);
            
            ydual = Q*results.x(1:K+2);
            solverTime = results.runtime;
            %fprintf('Qurobi QP time is : %16e \n', solverTime);
            
        otherwise
            disp('other value')
    end
end

%% local wrapper function for ASP
function [ydual, solverTime, inform] = ASPwrap(A,b,Q,tau,K)
    % ASP
    % min_y 1/2 lambda y'y -b'y 
    % subject to bl <= A'y <= bu
    
    ASPtimeVal = tic();
    [m,n] = size(A);
    
    % formulation transformation
    delta = 1e-6;

    D = [eye(K+2) zeros(K+2,1); zeros(1,K+2) delta];
    quadb = [Q'*b; -tau];
    Dinv = [eye(K+2) zeros(K+2,1); zeros(1,K+2) 1/delta];
    tildeb = Dinv*quadb; %since it's diagonal, D^(-1)=1/D  
    quadAprime = [A'*Q -ones(n,1); -A'*Q -ones(n,1)];
    tildeA = Dinv*quadAprime';  
    
    bl = zeros(2*n,1);
    bl(:) = -inf;
    bu = zeros(2*n,1);
    lambdain = 1;
    
    options = as_setparms;
    options.loglevel = 0;
    options.gaptol = 1e-12;
    [active,state,x,tildev,S,R,inform] = BPdual(tildeA,tildeb,bl,bu,lambdain,[], [], [], [], [], options);
    tildev
    %fprintf('hhhhhhh %15e \n ',tildev);
    % inverse transformation 
    v = Dinv * tildev; 
    ydual = v(1:K+2); 
    solverTime = toc(ASPtimeVal);
end