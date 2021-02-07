%% test rp_iso.m and rp_perp.m

disp('TEST rp_iso.m and rp_perp.m')

addpath('../.')

disp('- load data and test results')
x = load('data/noisesmooth_embed.dat');
rp_test_iso = load('data/noisesmooth_rpiso.dat');
rp_test_iso4 = load('data/noisesmooth_rpiso.4.dat');
rp_test_perp = load('data/noisesmooth_rpperp.dat');
rp_test_perp3 = load('data/noisesmooth_rpperp.3.dat');


disp('- rp_iso test 1')
r = rp_iso(x,.2,.3);
assert(isequal(rp_test_iso,r))
disp('  > passed')

disp('- rp_iso test 2')
r = rp_iso(x,.4,.2);
assert(isequal(rp_test_iso4,r))
disp('  > passed')

disp('- rp_perp test 1')
r = rp_perp(x,.2,.5);
assert(isequal(rp_test_perp,r))
disp('  > passed')

disp('- rp_perp test 2')
r = rp_perp(x,.3,.3);
assert(isequal(rp_test_perp3,r))
disp('  > passed')

disp('TEST all tests passed')
