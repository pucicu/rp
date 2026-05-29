%% test dl_censi.m

addpath('../.')

disp('- load test data')

x1 = load('data/noisesmooth_embed.dat');
x2 = embed(sin(linspace(0,2*pi*10,1000)),2,25);

dl_test_noise = load('data/noise_dl_censi.dat');
dl_test_sine = load('data/sine_dl_censi.dat');



disp('- dl_censi(noise)')

r = rp(x1,.1,'fix');

[dl, set_lines] = dl_censi(r);

y = [dl(:); set_lines(:)];

assert(isequal(round(100000*dl_test_noise), round(100000*y)))

disp('  > passed')



disp('- dl_censi(sine)')

r = rp(x2,.1,'fix','euc','vector');

[dl, set_lines] = dl_censi(r);

y = [dl(:); set_lines(:)];

assert(isequal(round(100000*dl_test_sine), round(100000*y)))

disp('  > passed')


disp('TEST all tests passed')
