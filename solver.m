function [ydual, solverTime] = solver(mode, A, b, p, Q, tau, varargin)
%-------------------------------------------------------------
% Provide different solvers for user to solve dual problem
% mode 1: built-in quadprog
% mode 2: pdco by Michael Saunder
% mode 3: Gurobi - commercial optimization toolbox
% mode 4: dnopt by Philip E Gill, Michael A Saunders and Elizabeth Wong
%-------------------------------------------------------------

    [m,n] = size(A);
    switch mode
        case 1
            % quadprog with "interior-point-convex" algorithm
            % min_x 1/2 xHx+c'x 
            % subject to Ax<=b or Aeq x = beq or lb <= x <= up
            
            cvec = [Q'*b;-tau];
            C = [A'*Q -ones(n,1);-A'*Q -ones(n,1)]; 
            bu = zeros(2*n,1); 
            
            quadTimeVal = tic();
            [x,fquad,exitflag] = quadprog(varargin{1},cvec,C,varargin{2},...
                                          [],[],[],[],[],varargin{3}); 
            solverTime = toc(quadTimeVal);
            ydual = Q*x(1:p+2,1);
          
            
        case 2
            % pdco (Primal-Dual Barrier Method)
            %    minimize    phi(x) + 1/2 norm(D1*x)^2 + 1/2 norm(r)^2
            %      x,r
            %    subject to  A*x + D2*r = beq,   bl <= x <= bu,   r unconstrained
            
            pdcoTimeVal = tic();
            cvec = [-Q'*b; tau; sparse(2*n,1)];
            C = [-A'*Q, -ones(n,1), speye(n), sparse(n,n); ...
                  A'*Q, -ones(n,1),sparse(n,n), speye(n)];
         
            [v,y,z,inform,PDitns,CGitns,time] = pdco(cvec,C,varargin{1},...
                varargin{2},varargin{3},varargin{5},1e-4,varargin{6},...
                varargin{4},sparse(2*n,1),sparse(p+2+1+2*n,1),1,1);
            solverTime   = toc(pdcoTimeVal);
            ydual = Q*v(1:p+2); 
            
        
        case 3
            % Gurobi
            GurobiTime = tic();
            
            Q1 = sparse(2*n+p+2+1,2*n+p+2+1);
            Q1(1:p+2,1:p+2) = 1/2.*Q'*Q;

            model.Q = Q1;
            model.obj = [-Q'*b; tau; zeros(2*n,1)];
            model.A = sparse([-A'*Q -ones(n,1) eye(n) zeros(n); ...
                               A'*Q -ones(n,1) zeros(n) eye(n)]);
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

        case 4 
            % dnopt
            dqoptTime = tic();
            
            cvec = [-Q'*b; tau];
            C = [-A'*Q -ones(n,1) ; A'*Q -ones(n,1)];
            
            [x,fval,exitFlag,output,lambda,states] = dqopt(varargin{1}, ...
                cvec,varargin{2}, varargin{3}, varargin{4}, C, ...
                varargin{5}, varargin{6},varargin{7});
            f = @(x) (cvec'*x + 1/2*x'*varargin{1}*x);
            g = @(x) varargin{1}*x+cvec;

            ydual = Q*x(1:p+2);
            solverTime = toc(dqoptTime);
            
        otherwise
            disp('other value')
    end
end

