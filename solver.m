function [ydual, solverTime] = solver(mode, A, b, p, Q, tau, arg)
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
            % min_x 1/2 xHx+c'x 
            % subject to Ax<=b or Aeq x = beq or lb <= x <= up
            
            cvec = [Q'*b;-tau];
            C = [A'*Q -ones(n,1);-A'*Q -ones(n,1)]; 
            bu = zeros(2*n,1); 
            
            quadTimeVal = tic();
            [x,fquad,exitflag] = quadprog(arg{1},cvec,C,arg{2},[],[],[],[],[],arg{3}); 
            solverTime = toc(quadTimeVal);
            ydual = Q*x(1:p+2,1);
          
            
        case 2
            % pdco (Primal-Dual Barrier Method)
            %    minimize    phi(x) + 1/2 norm(D1*x)^2 + 1/2 norm(r)^2
            %      x,r
            %    subject to  A*x + D2*r = beq,   bl <= x <= bu,   r unconstrained
            
            pdcoTimeVal = tic();
            %H = sparse(p+2+1+2*n, p+2+1+2*n); 
            %H(1:p+2, 1:p+2) = speye(p+2); %= Q'*Q is an identity matrix
            %objectivefunction = @(x) deal(0.5*(x'*H*x) + c'*x, H*x + c, H);
            cvec = [-Q'*b; tau; sparse(2*n,1)];
            C = [-A'*Q, -ones(n,1), speye(n), sparse(n,n); A'*Q, -ones(n,1),sparse(n,n), speye(n)];
           
            %include quadratic part into the obj function
            %[v,y,z,inform,PDitns,CGitns,time]=pdco(objectivefunction,Amatrix,bvec,bl,bu,1e-4,1e-4,options1,v0,sparse(2*n,1),sparse(p+1+2*n,1),1,1);
            %include quadratic part into the regularization D1
            [v,y,z,inform,PDitns,CGitns,time] = pdco(cvec,C,arg{1},arg{2},arg{3},...
                arg{5},1e-4,arg{6},arg{4},sparse(2*n,1),sparse(p+2+1+2*n,1),1,1);
            solverTime   = toc(pdcoTimeVal);
            ydual = Q*v(1:p+2); %v=[c,lambda, s1, s2], c: (p+2)x 1
            
            
        case 3
            % ASP (active-set pursuit)
            % ---------------------------------
            % min_y 1/2 lambda y'y -b'y 
            % subject to bl <= A'y <= bu
            % ---------------------------------
            
            [ydual, solverTime, inform] = ASPwrap(A,b,Q,tau,p);
            ydual = Q*ydual;
            fprintf('ASP time is : %16e \n', solverTime);
            
        
        case 4
            % Gurobi
            GurobiTime = tic();
            
            Q1 = sparse(2*n+p+2+1,2*n+p+2+1);
            Q1(1:p+2,1:p+2) = 1/2.*Q'*Q;

            model.Q = Q1;
            model.obj = [-Q'*b; tau; zeros(2*n,1)];
            model.A = sparse([-A'*Q -ones(n,1) eye(n) zeros(n); A'*Q -ones(n,1) zeros(n) eye(n)]);
            model.rhs = zeros(2*n,1);
            model.sense = '=';
            lb = zeros(2*n+p+2+1,1); 
            lb(1:p+2+1) = -inf; 
            model.lb = lb;
            gurobi_write(model, 'qp.lp');
            params.outputflag = 0; 
            results = gurobi(model,params);
            
            ydual = Q*results.x(1:p+2);
            solverTime = toc(GurobiTime);
            
        case 5
            % qpopt (not recommanded)
            qpoptTime = tic();
            cvec = [-Q'*b; tau];
            C = [-A'*Q -ones(n,1) ; A'*Q -ones(n,1)];
            [x,obj,lambda,istate,iter,inform] = qpopt(C,cvec,arg{1},arg{2},...
                arg{3},arg{4},1);
            ydual = Q*x(1:p+2);
            solverTime = toc(qpoptTime);
            
        case 6 
            % dqopt
            dqoptTime = tic();
            
            cvec = [-Q'*b; tau];
            C = [-A'*Q -ones(n,1) ; A'*Q -ones(n,1)];
            
            options.printfile = '';
            options.screen = 'off';
            
            [x,fval,exitFlag,output,lambda,states] = dqopt(arg{1}, cvec,...
                arg{2}, arg{3}, arg{4}, C, arg{5}, arg{6},options);
            ydual = Q*x(1:p+2);
            solverTime = toc(dqoptTime);
            
        case 7 
            % ADMM
            
        otherwise
            disp('other value')
    end
end

%% local wrapper function for ASP
function [ydual, solverTime, inform] = ASPwrap(A,b,Q,tau,p)
    % ASP
    % min_y 1/2 lambda y'y -b'y 
    % subject to bl <= A'y <= bu
    
    ASPtimeVal = tic();
    [m,n] = size(A);
    
    % formulation transformation
    delta = 1e-6;

    D = [eye(p+2) zeros(p+2,1); zeros(1,p+2) delta];
    quadb = [Q'*b; -tau];
    Dinv = [eye(p+2) zeros(p+2,1); zeros(1,p+2) 1/delta];
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
    
    
    % inverse transformation 
    v = Dinv * tildev; 
    ydual = v(1:p+2); 
    solverTime = toc(ASPtimeVal);
end