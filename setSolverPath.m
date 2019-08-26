% Add DNOPT/DQOPT matlab directories to Matlab path.

function []= setSolverPath(solverNum)
   
    switch solverNum
        case 2
            % pdco_dir = /path/to/pdco/solver
            pdco_dir = '/Users/casiebao/Dropbox/Research/spgl1/pdco-master';
            addpath([pdco_dir,'/code'                 ],'-end');
            
        case 4
            % dnopt_dir = /path/to/dnopt/solver
            dnopt_dir = '/Users/casiebao/Dropbox/Research/spgl1/dnopt';
            addpath([dnopt_dir,'/matlab'              ],'-end');
            addpath([dnopt_dir,'/matlab/util'         ],'-end');
            addpath([dnopt_dir,'/matlab/precompiled'  ],'-end');
    end
end
