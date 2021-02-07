%% test rp.m

disp('TEST rp.m')

addpath('../.')

disp('- load data and test results')
x = load('data/noisesmooth_embed.dat');
rp_test_max = load('data/noisesmooth_rp_max.dat');
rp_test_euc = load('data/noisesmooth_rp_euc_fix.dat');
rp_test_fan = load('data/noisesmooth_rp_euc_fan.dat');
rp_test_var = load('data/noisesmooth_rp_euc_var.dat');


disp('- rp(fix,max,loops)')
r = rp(x,.1,'fix','max','loops');
assert(isequal(rp_test_max,r))
disp('  > passed')

disp('- rp(fix,max,vector)')
r = rp(x,.1,'fix','max','vector');
assert(isequal(rp_test_max,r))
disp('  > passed')

disp('- rp(fix,max,matlabvector)')
r = rp(x,.1,'fix','max','matlabvector');
assert(isequal(rp_test_max,r))
disp('  > passed')

disp('- rp(fix,euc,loops)')
r = rp(x,.1,'fix','euc','loops');
assert(isequal(rp_test_euc,r))
disp('  > passed')

disp('- rp(fan,euc,loops)')
r = rp(x,.05,'fan','euc','loops');
assert(isequal(rp_test_fan,r))
disp('  > passed')

disp('- rp(var,euc,loops)')
r = rp(x,.05,'var','euc','loops');
assert(isequal(rp_test_var,r))
disp('  > passed')

disp('TEST all tests passed')


