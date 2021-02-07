%% test embed.m

disp('TEST embed.m')

addpath('../.')

disp('- set parameters')
m = 3;
tau = 2;

disp('- set input vector')
x = 1:10;

disp('- set desired result')
y_test = [
     1     3     5;
     2     4     6;
     3     5     7;
     4     6     8;
     5     7     9;
     6     8    10;
];

disp('- apply embed.m')
y = embed(x,m,tau);

disp('- compare result')
assert(isequal(y, y_test));

disp('TEST passed')
