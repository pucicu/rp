%% test rp.m

x = load('tests/data/noise_embed.dat');
rp_test_max = load('tests/data/noise_rp_max.dat');
rp_test_euc = load('tests/data/noise_rp_euc_fix.dat');
rp_test_fan = load('tests/data/noise_rp_euc_fan.dat');
rp_test_var = load('tests/data/noise_rp_euc_var.dat');

r = rp(x,.1,'fix','max','loops');
assert(isequal(rp_test_max,r))

r = rp(x,.1,'fix','max','vector');
assert(isequal(rp_test_max,r))

r = rp(x,.1,'fix','max','matlabvector');
assert(isequal(rp_test_max,r))

r = rp(x,.1,'fix','euc','loops');
assert(isequal(rp_test_euc,r))

r = rp(x,.05,'fan','euc','loops');
assert(isequal(rp_test_fan,r))

r = rp(x,.05,'var','euc','loops');
assert(isequal(rp_test_var,r))


