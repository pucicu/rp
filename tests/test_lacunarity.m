%% test lacunarity.m

addpath('../.')

disp('- load test data')

x1 = load('data/noisesmooth_embed.dat');
x2 = sin(linspace(0,2*pi*10,1000))' > .5;

boxSize = 1:50;

lac_test_noise = load('data/noise_lacunarity.dat');
lac_test_sine = load('data/sine_lacunarity.dat');


disp('- lacunarity(noise)')

r = rp(x1,.1,'fix');

y = lacunarity(double(r), boxSize);

assert(isequal(round(100000*lac_test_noise), round(100000*y)))

disp('  > passed')



disp('- lacunarity(sine)')

r = rp(x2,.1,'fix','euc','vector');

y = lacunarity(double(r), boxSize);

assert(isequal(round(100000*lac_test_sine), round(100000*y)))

disp('  > passed')


disp('TEST all tests passed')
