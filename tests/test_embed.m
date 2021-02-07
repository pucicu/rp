%% test embed.m

m = 3;
tau = 2;
x = 1:10;

y_test = [
     1     3     5;
     2     4     6;
     3     5     7;
     4     6     8;
     5     7     9;
     6     8    10;
];

y = embed(x,m,tau);

assert(isequal(y, y_test))
