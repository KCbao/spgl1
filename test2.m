s = RandStream('mt19937ar','Seed',0);

A = randn(s,50,256);
b = randn(s,50,1);

options = [];
options.iterations = 1000;
[x,r,g,info] = spgl1(A,b,3,[],[],options);
