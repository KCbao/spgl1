%{
s = RandStream('mt19937ar','Seed',0);

A = randn(s,50,256);
b = randn(s,50,1);

options = [];
options.iterations = 1000;
[x,r,g,info] = spgl1(A,b,3,[],[],options);
%}
%
function spgdemo(interactive)

if nargin < 1 || isempty(interactive), interactive = true; end
    
    % Initialize random number generators 
    rand('state',0);
    randn('state',0);
    
    % -----------------------------------------------------------
    % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
    %
    %    minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...
    %
    % A: 50 x 128, x: sparse support k=14, m/n=40%
    %
    % -----------------------------------------------------------
    
    % Create random m-by-n encoding matrix and sparse vector
    m = 50; n = 128; 
    k = 14;
    [A,Rtmp] = qr(randn(n,m),0);
    A  = A';
    p  = randperm(n); p = p(1:k);
    x0 = zeros(n,1); x0(p) = randn(k,1);
 
  
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the underdetermined LASSO problem for   \n');
    fprintf('%%                                               \n');
    fprintf('%%   minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...\n');
    fprintf('%%                                               \n');
    fprintf('%%  A: 50 x 128 sparse support k=14, m/n=0.4 \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);
    

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end

    % Set up vector b, and run solver    
    b = A * x0;
    tau = pi;
    options = [];
    options.iterations = 1000;
    [x,r,g,info] = spgl1(A,b,tau,[],[],options);

    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf(['nonzeros(x) = %i,   ' ...
             '||x||_1 = %12.6e,   ' ...
             '||x||_1 - pi = %13.6e\n'], ...
            length(find(abs(x)>1e-5)), norm(x,1), norm(x,1)-pi);
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');
    
    % -----------------------------------------------------------
    % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
    %
    %    minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...
    %
    % A: 50 x 256, x: sparse support k=14 m/n=19.5%
    %
    % -----------------------------------------------------------
    
    % Create random m-by-n encoding matrix and sparse vector
    
    m = 50; n = 256; 
    k = 14;
    [A,Rtmp] = qr(randn(n,m),0);
    A  = A';
    p  = randperm(n); p = p(1:k);
    x0 = zeros(n,1); x0(p) = randn(k,1);
 
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the underdetermined LASSO problem for   \n');
    fprintf('%%                                               \n');
    fprintf('%%   minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...\n');
    fprintf('%%                                               \n');
    fprintf('%%  A: 50 x 256 sparse support k=14, m/n=0.195 \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);
    

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end

    % Set up vector b, and run solver    
    b = A * x0;
    tau = pi;
 
    [x,r,g,info] = spgl1(A,b,tau,[],[],options);

    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf(['nonzeros(x) = %i,   ' ...
             '||x||_1 = %12.6e,   ' ...
             '||x||_1 - pi = %13.6e\n'], ...
            length(find(abs(x)>1e-5)), norm(x,1), norm(x,1)-pi);
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');
    
    % -----------------------------------------------------------
    % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
    %
    %    minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...
    %
    % A: 50 x 512, x: sparse support k=14 m/n=1%
    %
    % -----------------------------------------------------------
    
    % Create random m-by-n encoding matrix and sparse vector
    
    m = 50; n = 512; 
    k = 14;
    [A,Rtmp] = qr(randn(n,m),0);
    A  = A';
    p  = randperm(n); p = p(1:k);
    x0 = zeros(n,1); x0(p) = randn(k,1);
 
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the underdetermined LASSO problem for   \n');
    fprintf('%%                                               \n');
    fprintf('%%   minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...\n');
    fprintf('%%                                               \n');
    fprintf('%%  A: 50 x 512 sparse support k=14, m/n=0.098 \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);
    

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end

    % Set up vector b, and run solver    
    b = A * x0;
    tau = pi;
    
    [x,r,g,info] = spgl1(A,b,tau,[],[],options);

    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf(['nonzeros(x) = %i,   ' ...
             '||x||_1 = %12.6e,   ' ...
             '||x||_1 - pi = %13.6e\n'], ...
            length(find(abs(x)>1e-5)), norm(x,1), norm(x,1)-pi);
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');
    
    % -----------------------------------------------------------
    % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
    %
    %    minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...
    %
    % A: 50 x 512, x: sparse support k=20 m/n=1%
    %
    % -----------------------------------------------------------
      % Create random m-by-n encoding matrix and sparse vector
    
    m = 50; n = 512; 
    k = 20;
    [A,Rtmp] = qr(randn(n,m),0);
    A  = A';
    p  = randperm(n); p = p(1:k);
    x0 = zeros(n,1); x0(p) = randn(k,1);
 
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the underdetermined LASSO problem for   \n');
    fprintf('%%                                               \n');
    fprintf('%%   minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...\n');
    fprintf('%%                                               \n');
    fprintf('%%  A: 50 x 512 sparse support k=20 , m/n=0.098 \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);
    

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end

    % Set up vector b, and run solver    
    b = A * x0;
    tau = pi;
    
    [x,r,g,info] = spgl1(A,b,tau,[],[],options);

    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf(['nonzeros(x) = %i,   ' ...
             '||x||_1 = %12.6e,   ' ...
             '||x||_1 - pi = %13.6e\n'], ...
            length(find(abs(x)>1e-5)), norm(x,1), norm(x,1)-pi);
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');
    
     % -----------------------------------------------------------
    % Solve the underdetermined LASSO problem for ||x||_1 <= pi:
    %
    %    minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...
    %
    % A: 50 x 512, m/n=1%, b completely random
    %
    % -----------------------------------------------------------
      % Create random m-by-n encoding matrix and sparse vector
    
    m = 50; n = 512; 
    k = n;
    [A,Rtmp] = qr(randn(n,m),0);
    A  = A';
    p  = randperm(n); p = p(1:k);
    x0 = zeros(n,1); x0(p) = randn(k,1);
 
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the underdetermined LASSO problem for   \n');
    fprintf('%%                                               \n');
    fprintf('%%   minimize ||Ax-b||_2 subject to ||x||_1 <= 3.14159...\n');
    fprintf('%%                                               \n');
    fprintf('%%  A: 50 x 512, m/n=0.098, b completely random \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);
    

    fprintf('\nPress <return> to continue ... \n');
    if interactive, pause; end

    % Set up vector b, and run solver    
    b = A * x0;
    b
    randn(50,1)
    tau = pi;
    
    [x,r,g,info] = spgl1(A,b,tau,[],[],options);

    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf(['nonzeros(x) = %i,   ' ...
             '||x||_1 = %12.6e,   ' ...
             '||x||_1 - pi = %13.6e\n'], ...
            length(find(abs(x)>1e-5)), norm(x,1), norm(x,1)-pi);
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');
    
    
end % function demo 

%}