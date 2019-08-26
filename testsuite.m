rng('default') % set a pointer to 0
% Choose a test problem by problem index
index = 3; 
[A, b, tau, options] = generateProblem(index);
options.dualMode = 1;
% options.hybridMode = 1; 
% original spgl1
x = spgl1(A,b,tau,[],[],options); 
%% Private function generateProblem
function [A, b, tau, options]=generateProblem(index)
  
  options = [];
  options.iterations = 1e+6;
      
  switch index
    case 1
        % -----------------------------------------------------------
        % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
        %
        %    minimize ||Ax-b||_2 subject to ||x||_1 <= tau...
        %
        % A: 50 x 128
        %
        % -----------------------------------------------------------
       
        m = 50; n = 256;
        A = randn(m,n);
        s = 10; % sparsity
        x = zeros(n,1);
        p = randperm(n,s);
        x(p) = sign(randn);
        b = A*x;
        tau= 0.1*norm(x,1);

        fprintf(['%% ', repmat('-',1,78), '\n']);
        fprintf('%% Solve the underdetermined LASSO problem for   \n');
        fprintf('%%                                               \n');
        fprintf('%%   minimize ||Ax-b||_2 subject to ||x||_1 <= ~ taubp...\n');
        fprintf('%%                                               \n');
        fprintf('%%  A: 50 x 128 with no sparse support on x \n');
        fprintf(['%% ', repmat('-',1,78), '\n']);

        % original solver 2695 iterations
        % our version: 778 iters
        % final average mutual coherence: 0.336477
        % just append Rr: 965 iters, mutual coh: 0.498796

      case 2
        % -----------------------------------------------------------
        % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
        %
        %    minimize ||Ax-b||_2 subject to ||x||_1 <= tau...
        %
        %  A: 50 x 256
        %
        % -----------------------------------------------------------
        s = RandStream('mt19937ar','Seed',0);
        m = 50; n = 256;
        A = randn(s,m,n);
        b = randn(s,m,1);

        tau=4.2071091e+00;
        %4.2085938e+00 with sigma=1e-3 3201 iters
        %spg_bp(A, b, opts);

        fprintf(['%% ', repmat('-',1,78), '\n']);
        fprintf('%% Solve the underdetermined LASSO problem for   \n');
        fprintf('%%                                               \n');
        fprintf('%%   minimize ||Ax-b||_2 subject to ||x||_1 <=  ~ taubp...\n');
        fprintf('%%                                               \n');
        fprintf('%%  A: 50 x 256 with no sparse support on x \n');
        fprintf(['%% ', repmat('-',1,78), '\n']);

        % original solver 3454 iterations
        % set tau=3, takes 145 iterations
        % set tau=4.5 > tau_{BP}, takes 115 iterations
        % when very close to tau_{BP}, it takes the most num of iters
        % our version: 931 iters with stopping criteria: rgap_accel<optTol only
        % final average mutual coherence: 0.349046
        % just append Rr: 1097 iters, mutual coh:  0.499660





    
      case 3
        % -----------------------------------------------------------
        % Solve
        %
        %    minimize ||y||_1 subject to AW^{-1}y = b
        %    
        % and the weighted basis pursuit (BP) problem:
        %
        %    minimize ||Wx||_1 subject to Ax = b
        %
        % followed by setting y = Wx.
        % -----------------------------------------------------------

        % Initialize random number generators 
        rand('state',3);
        randn('state',7);

        % Create random m-by-n encoding matrix and sparse vector
        m = 50; n = 128; k = 14;
        [A,Rtmp] = qr(randn(n,m),0);
        A  = A';
        p  = randperm(n); 

        % Sparsify vector x0 a bit more to get exact recovery
        k = 9;
        x0 = zeros(n,1); x0(p(1:k)) = randn(k,1);

        % Set up weights w and vector b 
        w     = rand(n,1) + 0.1;       % Weights
        b     = A * (x0 ./ w);         % Signal

        % Run solver for both variants
        
        AW   = A * spdiags(1./w,0,n,n);
        %tau = 15;
        tau = 1.04010450e+01; %hard-code ones, original 1.0401468e+01

        A=AW;

        fprintf(['%% ', repmat('-',1,78), '\n']);
        fprintf('%% Solve                                          \n');
        fprintf('%%                                                \n');
        fprintf('%% (1) minimize  AW^{-1}y = b  subject to ||y||_1 <= taubp \n');
        fprintf('%%                                                \n');
        fprintf(['%% ', repmat('-',1,78), '\n']);

        % original solver: 1407 iters
        % our solver pdco: 636 iters  7.0943022e-05 
        % total time: 15.9s pdco: 14.3s
        % average mutual coherence 0.318467
        % just append Rr: 638 iters, mutual coh: 0.488968
        
        % Gurobi: 598 iters  9.1668574e-05 total: 10.9s Gurobi: 3.4s
        % ADMM: 598 iters, 9.1668574e-05, total 5.4s, 4.9s and 500 iters in
        % ADMM
        % qpopt : 598 iters, 9.1668574e-05 total: 0.9s, solver: 0.2s
        
        % evrey 20: 9.3658740e-05, 600, 0.7s, 0.2s
        % Gurobi 10: 605 1.2 1.0 8.6894617e-05 
        % Gurobi 20: 615 0.7 0.5 8.4234041e-05
        % Gurobi with 40: 615 0.5 0.2  8.4234041e-05 
        % Gurobi with (1+alpha)^k:  627 1.5 1.3 7.7803209e-05
    
    case 4
        % Iterations with opttol 1e-4
        % v2 Ewout Original: 47366 (148.1 sec.) Hybrid  : 47366 (151.4 sec.)
        % v2 mine original: 45437 36.7s, hybrid: 45437 37.5s
        
        % Gurobi 10: 45665 867.9 827.5 2.7478385e-04
        % Gurobi 20: 45795 478.7 438.3 2.7477840e-04
        % Gurobi with 40: 45815 252.5 218.2
        % Gurobi with (1+alpha)^k:  46941 (73.1 29.6 2.7473821e-04)
        s = RandStream('mt19937ar','Seed',1);
        A = randn(s,512,1024);
        b = randn(s,512,1);
        tau = 23; % Sigma = 2.344e-2


    case 5 % (case 1 with smaller tau)
        % -----------------------------------------------------------
        % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
        %
        %    minimize ||Ax-b||_2 subject to ||x||_1 <= tau...
        %
        %  A: 512 x 1024, b: 512 x 1, tau = 22
        %  opttol: 1e-4
        % -----------------------------------------------------------
      
        % v1 Mine original: 9672 8.6s 
        % v2 Ewout Original: 8757 (27.9 sec.) Hybrid  : 1805 ( 6.0 sec.)
        % v2 Mine original: 10079 8.3s, hybrid: 2010 2.2s
        
        
        % pdco: 3590 (total:  829.9s, pdco:  806.6s),average mutual coherence is u = 0.500000
        % Gurobi: 3191 (total: 579.3s, Gurobi: 168.5s), average mutual coherence is u = 0.500000
        % Gurobi sol: 2.4417391e-01
        
        % dpopt: 4640: 76.9s, 69.6s 2.4417365e-01
        % with warm-start: 4557: 50s 43.1s
        
        % Gurobi 10: 3295 63.7 60.6 2.4417389e-01
        % Gurobi 20: 4595 46.1 42.4 2.4417366e-01
        % Gurobi with 40: 4135 22.0 18.8
        % Gurobi with (1+alpha)^k: 5911 26.6 22.1 2.4417351e-01 
        s = RandStream('mt19937ar','Seed',1);
        A = randn(s,512,1024);
        b = randn(s,512,1);
        tau = 22; % Sigma = 6.988e-1

    % ---

    case 6
        % -----------------------------------------------------------
        % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
        %
        %    minimize ||Ax-b||_2 subject to ||x||_1 <= tau...
        %
        %  A: 256 x 1024, b: 256 x 1, tau = 11.1
        %  opttol: 1e-4
        % -----------------------------------------------------------
 
        % v1 Mine: 42409 20.7s 1.3001443e-02
        % v2 Ewout Original: 40332 (56.4 sec.)  Hybrid  : 40332 (59.4 sec.) 
        % v2 Mine original: 39822 16.8s, hybrid: 39822 17.4s
        
        % 
        % pdco:  40829  8.4682352e-05 
        % Gurobi: 15051 (2958.6s, 838.1s)average mutual coherence is u = 0.418855
        % sol:  9.9999022e-05
        
        % dqopt: 15051 (157.9, 141.1)
        % with warm-start: 15051 131.6 115.0
        % sol: 9.9999022e-05
        
        % Gurobi 10: 15055 270.9 259.8 9.9991875e-05
        % Gurobi 20: 15055 140.1 131.1
        % Gurobi with 40: 15055 72.6 64.8 9.9991875e-05
        % Gurobi with (1+alpha)^k:  15685 33.2 25.1 9.8764713e-05
        s = RandStream('mt19937ar','Seed',3);
        A = randn(s,256,1024);
        b = randn(s,256,1);
        tau = 11.1;

    case 7 % (case 3 with smaller tau)
        % -----------------------------------------------------------
        % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
        %
        %    minimize ||Ax-b||_2 subject to ||x||_1 <= tau...
        %
        %  A: 256 x 1024, b: 512 x 1, tau = 10.9
        %  opttol: 1e-4
        % -----------------------------------------------------------
     
        % v1 Mine original: 5599 3.0s) 2.5027851e-01
        % v2 Ewout Original: 5175 (7.2 sec) Hybrid  : 3169 (5.1 sec.)
        % v2 Mine Original: 4441 (2.2s)  Hybrid  : 3740 (2.2s)
        
        % pdco: 3872 (total time:  889.5s, pdco cost: 873.4s)
        % Gurobi: 3642  3.1320725e-02, total:658.1s, gurobi:
        % 185.1s
        
        % dqopt: 4274: 51.1s, 46.3s
        % with warm-start: 4376 43.1s 36.3s
        
        % Gurobi 10: 3775 69 66.1 3.1320437e-02
        % Gurobi 20: 3935 37.7 35.2 3.1320292e-02
        % Gurobi with 40: 4095 20.6 18.3  3.1320138e-02
        % Gurobi with (1+alpha)^k: 4632, 22.8, 20.4 3.1319853e-02 
        
        s = RandStream('mt19937ar','Seed',3);
        A = randn(s,256,1024);
        b = randn(s,256,1);
        tau = 10.9;


    case 8 % (case 3 with smaller tau)
        % -----------------------------------------------------------
        % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
        %
        %    minimize ||Ax-b||_2 subject to ||x||_1 <= tau...
        %
        %  A: 256 x 1024, b: 256 x 1, tau = 10.7
        %  opttol: 1e-4
        % -----------------------------------------------------------
        
        % Iterations with opttol 1e-4
        % v1 Mine original: 2168 (1.4s)
        % v2 Ewout Original: 2028 (2.7 sec.) Hybrid  : 1699 (2.7 sec.)
        % v2 Mine original: 2168 (1.1s) Hybrid: 1723 (1.1s)
        
        % pdco: 1987 (total time: 466.9s, pdco: 458.6s)
        % solution norm(x,1) = 10.7
        % Gurobi: 1970 (total time: 351.6s, Gurobi: 101.7s)
        % dqopt : 1986 (total time: 24.8, solver: 22.5)
        % dqopt: with warm-start: 1981 18.4 16.1
        % norm(x,1)
        
        % Gurobi 10: 1995 36.9 35.4 1.1914999e-01
        % Gurobi 20: 1995, 19.3,  (1.1914999e-01)
        % Gurobi with 40: 1995 18.9 17.6 1.1914999e-01  
        % Gurobi with (1+alpha)^k: 2021, 18.8, 17.7,  1.1914999e-01
        
        s = RandStream('mt19937ar','Seed',3);
        A = randn(s,256,1024);
        b = randn(s,256,1);
        tau = 10.7;

    % ---

    case 9
        % -----------------------------------------------------------
        % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
        %
        %    minimize ||Ax-b||_2 subject to ||x||_1 <= tau...
        %
        %  A: 512 x 768, b: 512 x 1, tau = 28.8
        %  opttol: 1e-4
        % -----------------------------------------------------------
        
        % v1 Mine original: 
        % v2 Ewout Original: 64010 (127.7 sec.) Hybrid  : 64010 (131.4 sec.)
        % v2 Mine Original: 64399 (34.7s ) Hybrid: 64399 (33.6s)
    
        
        % 63848 p=8, 52.6s
        
        % Gurobi: 61150, 7516.6s, solver time: 4964.3s
        
        % dqopt: 63424, 717.1s, 634.5 2.1795864e-04 
        
        % Gurobi 10: 62165 788.8 745.3 2.1802616e-04
        % Gurobi with 20: 63335 423.1 382.9
        % Gurobi with 40: 63335 222.7 186.0 2.1795969e-04
        % Gurobi with (1+alpha)^k: 64560, 59.2, 19.7,  2.1793189e-04
        s = RandStream('mt19937ar','Seed',6);
        A = randn(s,512,768);
        b = randn(s,512,1);
        tau = 28.8;

    case 10
        % -----------------------------------------------------------
        % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
        %
        %    minimize ||Ax-b||_2 subject to ||x||_1 <= tau...
        %
        %  A: 512 x 768, b: 512 x 1, tau = 28.0
        %  opttol: 1e-4
        % -----------------------------------------------------------
        
        % v2 Ewout Original: 22572 (45.2 sec.) Hybrid  :  8585 (19.0 sec.) 
        % v2 Mine Original: 24086 (12.0s) Hybrid: 8710 (5.3s)
        
        % our with pdco: 14519 (total time: 2473.6s, pdco:  2402.8s
        % average mutual coherence is u = 0.499999
        
        % Gurobi: 14205(total: 1781.3s, Gurobi 608.2s)
        % u = 0.499999 6.4487515e-02 
        
        % qpopt: 15586 199.9s, 177.3s
        % with warm-start: 15971 141.7s 121.4s
        
        % Gurobi 10: 14345, 183.2, 173.2
        % Gurobi with 20: 14515,  95.4, 86.5 6.4487430e-02
        % Gurobi with 40: 15655 56 46.8, 6.4487162e-02
        % Gurobi with (1+alpha)^k: 20018, 29.7, 17.2,  6.4486639e-02
        s = RandStream('mt19937ar','Seed',6);
        A = randn(s,512,768);
        b = randn(s,512,1);
        tau = 28.0;

    case 11
        % -----------------------------------------------------------
        % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
        %
        %    minimize ||Ax-b||_2 subject to ||x||_1 <= tau...
        %
        %  A: 512 x 768, b: 512 x 1, tau = 27.4
        %  opttol: 1e-4
        % -----------------------------------------------------------
        
        % v2 Ewout Original: 6871 (13.8 sec.) Hybrid  : 3048 ( 7.0 sec.) 
        % v2 Mine Original: 5891 (3.5s) Hybrid: 3272 (2.3s)
        
        % pdco: 6021  (total time:  1068.0s, pdco: 1038.8 s
        % average mutual coherence is u = 0.500000
        
        % Gurobi qp: 5392 (total time: 668s, Gurobi: 652.2, solver Gurobi: 219s)
        % average mutual coherence is u = 0.49999
        % sol: 1.8909e-01
        
        % qpopt: 5998 71.4s 63.5s
        % with warm-start: 6018 52.6 44.8
        
        % Gurobi 10: 6025 76.2 72.0 1.8909454e-01
        % Gurobi 20: 6055 (39.9, 36.1 ) 1.8909452e-01
        % Gurobi 40: 6055 (21.7, 18.1) 1.8909452e-01 
        % Gurobi with (1+alpha)^k: 7185 (19, 14.8) 1.8909425e-01 
        s = RandStream('mt19937ar','Seed',6);
        A = randn(s,512,768);
        b = randn(s,512,1);
        tau = 27.4;

        
        
 end
end




