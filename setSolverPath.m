% Add DNOPT/DQOPT matlab directories to Matlab path.

function []= setSolverPath(solverNum, solver_path_dir)
   
    switch solverNum
        case 2
            % pdco_dir = /path/to/pdco/solver
%             pdco_dir = '/Users/casiebao/Dropbox/Research/spgl1/pdco-master';
            addpath([solver_path_dir,'/code'                 ],'-end');
            
        case 4
            % dnopt_dir = /path/to/dnopt/solver
%             dnopt_dir = '/Users/casiebao/Dropbox/Research/spgl1/dnopt';
            addpath([solver_path_dir,'/matlab'              ],'-end');
            addpath([solver_path_dir,'/matlab/util'         ],'-end');
            addpath([solver_path_dir,'/matlab/precompiled'  ],'-end');
    end
end
