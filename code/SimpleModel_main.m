% Master script to run other scripts in desired order
tic
warning('off', 'MATLAB:polyshape:repairedBySimplify');
warning('off', 'MATLAB:polyshape:boundary3Points');
run('SimpleModel_initialize.m'); % complete
run('SimpleModel_inversionmap.m'); % complete
run('SimpleModel_safe_targ_recur.m'); % complete
run('SimpleModel_EXtraj.m'); % complete
run('SimpleModel_reachability.m'); % complete
run('SimpleModel_switchregions.m'); % complete
run('SimpleModel_OPTtraj.m'); % complete
toc