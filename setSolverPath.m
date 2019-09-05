% Add DNOPT/DQOPT matlab directories to Matlab path.

function []= setSolverPath(solverNum, solver_path_dir)
   
    switch solverNum
        case 2
            addpath([solver_path_dir,'/code'                 ],'-end');
            
        case 4
            addpath([solver_path_dir,'/matlab'              ],'-end');
            addpath([solver_path_dir,'/matlab/util'         ],'-end');
            addpath([solver_path_dir,'/matlab/precompiled'  ],'-end');
    end
end
