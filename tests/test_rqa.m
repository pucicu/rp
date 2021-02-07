%% test rqa.m

addpath('../.')

disp('- load data and test results')
x1 = load('data/noisesmooth_embed.dat');
x2 = embed(sin(linspace(0,2*pi*10,1000)),2,25);
x3 = [-0.0855, 2.3733, -0.4739, 0.9463, 0.8182, 1.5890, 0.5260, -2.4652, -0.8525, 0.5117]';
rqa_test_noisesmooth = load('data/noisesmooth_rqa.dat');
rqa_test_sine = load('data/sine_rqa.dat');
rqa_test_noise = load('data/noise_rqa.dat');

disp('- rqa(smoothed noise,fix,max,vector)')
r = rp(x1,.1,'fix','max','vector');
y = rqa(r,5,3);
assert(isequal(round(1000*rqa_test_noisesmooth),round(1000*y)))
disp('  > passed')

disp('- rqa(sin,fix,euc,vector)')
r = rp(x2,.1,'fix','euc','vector');
y = rqa(r,2,1);
assert(isequal(round(1000*rqa_test_sine),round(1000*y)))
disp('  > passed')




disp('- rqa(noise,fan)')
r = rp(x3,.5,'fan');
y = rqa(r,2,2);
assert(isequal(round(1000*rqa_test_noise),round(1000*y)))
disp('  > passed')


disp('TEST all tests passed')
