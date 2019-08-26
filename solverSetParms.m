function argout = solverSetParms(mode, p, A)
% ----------------------
% This function is to set one-time setup for some arguments for solver
% 
% ----------------------
    
     [m,n] = size(A);
     switch mode
         case 1
             % quadprog setup
             H = eye(p+2+1,p+2+1);
             H(p+2+1,p+2+1) = 0;
             bu = zeros(2*n,1); 
             argout = cell(3,1);
             options = optimoptions(@quadprog,'Display','off');
             argout{1} = H;
             argout{2} = bu;
             argout{3} = options;
         
         case 2 
             % pdco setup    
             beq = sparse(2*n, 1); 
             bl = [-Inf(p+2+1,1);sparse(2*n,1)]; 
             bu = Inf(p+2+1+2*n,1);
             v0 = [0;zeros(p+2+2*n,1)]; % warm start
             d1 = zeros(p+2+1+2*n,1); 
             d1(1:p+2) = 1;
             options = pdcoSet;
             options.Method = 22;
             options.Print = 0; % no print

             argout = cell(6,1);
             argout{1} = beq;
             argout{2} = bl;
             argout{3} = bu;
             argout{4} = v0;
             argout{5} = d1;
             argout{6} = options;
             
             
         case 3
             % Gurobi
             argout{1} = [];
         case 4
             % dnopt setup
             H = eye(p+2+1,p+2+1);
             H(p+2+1,p+2+1) = 0;
    
             x0 = zeros(p+2+1,1); % initialization
             x0(p+2+1) = 1;
             xl = -inf(p+2+1,1);
             xu = inf(p+2+1,1);
             bu = zeros(2*n,1);
             bl = -inf(2*n,1);
             
             options.printfile = '';
             options.screen = 'off';
             options.optimality_tolerance = 1e-06;
             options.iterations_limit = 100;  
             
             argout{1} = H;
             argout{2} = x0;
             argout{3} = xl;
             argout{4} = xu;
             argout{5} = bl;
             argout{6} = bu;
             argout{7} = options;
             
         otherwise
            disp('other value')
    end
    

end