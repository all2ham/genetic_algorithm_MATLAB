% Usage: gasimple_multi or gasimple_multi('x*exp(-x)');

function [bestsol, bestfun, count]=gasimple_edit(funstr)
global solnew sol pop popnew fitness fitold f range dim nsbit;
if nargin<1,
    % defining new function
    dim = 30; %%%%%%%%%%%%%% set desired amount of dimentions here
    funstr = '0';
    for i = 1:dim
        newterm = ['+ (abs(x('  num2str(i)  ')))^(' num2str(i) '+1)'];
        funstr = strcat(funstr, newterm);
    end
    % Easom Function with fmax=l at x=pi
    %funstr='for i=1:n;(abs(x))^2';
end
range=[-1 1]; % Range/Domain
% Converting to an inline function
f=vectorize(inline(funstr));
% Initializing the parameters
rand('state' ,0'); % Reset the random generator
popsize=20; % Population size
MaxGen=100; % Max number of generations
count=0;    % counter
nsite=2;    % number of mutation sites
pc=0.95;    % Crossover probability
pm=0.05;    % Mutation probability
nsbit=16;   % String length (bits)
% Generating the initial population
if dim ==1
    popnew=init_gen(popsize,nsbit);
else
    popnew=init_gen_multivariable(popsize,nsbit,dim);
end

fitness=zeros(1,popsize); % fitness array

% Initialize solution <- initial population
solnew = zeros(popsize,dim);
for i=1:popsize,
    solnew(i,:)=bintodec_multi(popnew(i,:),dim,nsbit);
end
% Start the evolution loop
for i=1:MaxGen,
    % Record as the history
   fitold=fitness; 
   pop=popnew; 
   sol=solnew;
   for j=1:popsize,
       % Crossover pair
       ii=floor(popsize*rand)+1; 
       jj=floor(popsize*rand)+1;
       % Cross over
       if pc>rand,
          [popnew(ii,:),popnew(jj,:)]=crossover(pop(ii,:),pop(jj,:));
          % Evaluate the new pairs
           count=count+2;
           evolve(ii); 
           evolve(jj);
       end
       % Mutation at n sites
       if pm>rand,
              kk=floor(popsize*rand)+1; 
              count=count+1;
              popnew(kk,:)=mutate(pop(kk,:),nsite);
              evolve(kk);
       end
    end % end for j
    % Record the current best
    bestfun(i)=max(fitness);
    bestsol(i)=mean(sol(bestfun(i)==fitness));
end

% Display results
set(gcf,'color','w');
subplot (2,1,1); plot(bestsol); title('Best estimates');
subplot(2,1,2); plot(bestfun); title('Fitness');



%%%%%%%%%%%%%%%%% All the sub functions %%%%%%%%%%%%%%%%%


% generation of the initial population
function pop=init_gen(np,nsbit)
% String length=nsbit+l with pop(:,l) for the Sign
pop=rand(np,nsbit+1)>0.5;
function pop=init_gen_multivariable(np,nsbit,dim)
%modified form of init_gen which composes strings for multiple variables
if dim == 1
    warning('You are using the multivariable version of init_gen.')
end
pop=rand(np,(nsbit+1)*dim)>0.5;
% Evolving the new generation
function evolve(j)
    global solnew popnew fitness fitold pop sol f dim nsbit;
    solnew(j,:)=bintodec_multi(popnew(j,:),dim,nsbit);
    fitness(j)=f(solnew(j,:));
    if fitness(j)>fitold(j),
        pop(j,:)=popnew(j,:);
        sol(j)=solnew(j);
    end
    
% Convert a binary string into a decimal number
function [dec]=bintodec(bin)
    global range;
    % Length of the string without sign
    nn=length(bin)-1;
    num=bin(2:end); % get the binary
    % Sign=+1 if bin(l)=0; Sign=-l if bin(l)=l.
    Sign=1-2*bin(1);
    dec=0;
    % floating point/decimal place in a binary string
    dp=floor(log2(max(abs(range))));
    for i=1:nn,
        dec=dec+num(i)*2^(dp-i);
    end
    dec=dec*Sign;
    
 function [dec]=bintodec_multi(bin,dim,bits)
    % modified version of bintodec to handle multivariable bitstrings
    global range;
    % Length of the string without signs
    % nn=length(bin)-dim;
    % get the binaries
    for i = 1:dim
        dtemp(i,:) = bin((bits*(i-1)+i+1):(bits*i+i));
        signs(i) = 1-2.*bin(bits*(i-1)+i);
    end
    % Sign=+1 if bin(l)=0; Sign=-l if bin(l)=l.
    dec=zeros(1,dim);
    % floating point/decimal place in a binary string
    dp=floor(log2(max(abs(range))));
    for j = 1:length(dec)
        for i=1:bits,
            dec(j)=dec(j)+dtemp(j,i)*2^(dp-i);
        end
    end
    dec=dec.*signs;
    
% Crossover operator
function [c,d]=crossover(a,b)

nn=length(a)-1;
% generating a random crossover point
cpoint=floor(nn*rand)+1;
c=[a(1:cpoint) b(cpoint+1:end)];
d=[b(1:cpoint) a(cpoint+1:end)];
% Mutatation operator
function anew=mutate(a,nsite)
    nn=length(a); anew=a;
    for i=1:nsite,
    j=floor(rand*nn)+1;
    anew(j)=mod(a(j)+1,2);
    end
