function result = QPSO(fun,D,nPop,lb,ub,maxit)
% INPUT:
%   fun     : function handle for optimization
%   D       : problem dimension (number of variables)
%   nPop    : number of particles in the swarm
%   lb      : lower bound constrain
%   ub      : upper bound constrain
%   maxit   : max number of iterations
%   maxeval : max number of function evaluations
% OUTPUT:
%   xmin    : best solution found
%   fmin    : function value at the best solution, f(xmin)
%   histout : record of function evaluations and fitness value by iteration
% EXAMPLE:
fun = @Objective; %fonksiyon
nPop = 100;
lb = -8;
ub = 8;
maxit = 1000; %max_iteration
SingMul=1;   %% single objective =1, multiobjective =2
DtypeG=19;   %% 4;12;19
D= DtypeG*256; % her bir zaman serisindeki deðiþken sayýsý {1024,3072, 4868}  
coll=256;
%sonuc= zeros(1,maxit);
sonuc= zeros(30,1);
 %maxeval = 10000*D; 
%[xmin,fmin,histout] = QPSO(fun,D,nPop,lb,ub,maxit,maxeval);
% OR DIRECTLY:
% [xmin,fmin,histout] = QPSO(@griewankfcn,30,50,-600,600,1000,10000*30);
% QPSO parameters:
w1 = 0.5;
w2 = 1.0;
c1 = 1.5;
c2 = 1.5;
%%% datalar okundu %%%%%
Amatrix=dlmread('D19NA.txt');
ICAcomponent=dlmread('D19NS.txt');
Mixed=dlmread('D19NX.txt');

for mmm=1:30
% Initializing solution
x = unifrnd(lb,ub,[nPop,D]);
% Evaluate initial population
pbest = x;
histout = zeros(maxit,2);
f_x=fun(x,Amatrix,Mixed,ICAcomponent,DtypeG,SingMul,nPop,coll); %Calculate the fitness values of solutions
%f_x = fun(x,func_num);
fval = nPop;
f_pbest = f_x;
[~,g] = min(f_pbest);
gbest = pbest(g,:);
f_gbest = f_pbest(g);
it = 1;
%histout(it,1) = fval;
%histout(it,2) = f_gbest;
while it <= maxit 
    alpha = (w2 - w1) * (maxit - it)/maxit + w1;
    mbest = sum(pbest)/nPop;
    for i = 1:nPop
        fi = rand(1,D);
        p = (c1*fi.*pbest(i,:) + c2*(1-fi).*gbest)/(c1 + c2);
        u = rand(1,D);
        
        b = alpha*abs(x(i,:) - mbest);
        v = log(1./u);
        if rand < 0.5
            x(i,:) = p + b .* v;
        else
            x(i,:) = p - b .* v;
        end
        
        % Keeping bounds
        x(i,:) = max(x(i,:),lb);
        x(i,:) = min(x(i,:),ub);
       % f_x(i) = fun(x(i,:), func_num);
       fval = fval + 1;
    end
     f_x=fun(x,Amatrix,Mixed,ICAcomponent,DtypeG,SingMul,nPop,coll); %Calculate the fitness values of solutions

        for i = 1:nPop
        if f_x(i) < f_pbest(i)
            pbest(i,:) = x(i,:);
            f_pbest(i) = f_x(i);
        end
        if f_pbest(i) < f_gbest
            gbest = pbest(i,:);
            f_gbest = f_pbest(i);
        end
        end
    %sonuc(1,it)=f_gbest; 
    it = it + 1;
    %histout(it,1) = fval;
    %histout(it,2) = f_gbest;
end
xmin = gbest; 
fmin = f_gbest;
sonuc(mmm,1)=fmin;
end

result.sonucc = sonuc;
end