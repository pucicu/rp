%% test rp.m

x = load('data/noise_embed.dat');
rp_test_iso = load('data/noise_rpiso.dat');
rp_test_perp = load('data/noise_rpperp.dat');

r = rp_iso(x,.2,.3);
assert(isequal(rp_test_iso,r))

r = rp_perp(x,.2,.5);
assert(isequal(rp_test_perp,r))

