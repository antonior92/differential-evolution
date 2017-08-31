function [xBest, fBest, nGen, X, fx, gx ] = de(fitnessfcn, nVars, lowerBound, upperBound, constr, Options)  
% DE - Diferential Evolution Algorithm
% Antonio Horta Ribeiro - 2016
% Belo Horizonte - Brasil
%
% It is a global minimization algorithm proposed by Storn, R. and Price, K. [1].
%
% function [xBest, fBest, nGen, X, fx, gx ] = de(fitnessfcn, nVars, lowerBound, upperBound, constr, Options)   
% Looks for the value 'xbest' that minimize 'fitnessfct' subject to
% 'constr'.
%
% Inputs:
% -> 'fitnessfcn': Handle to the fitness function. The fitness
% function should accept a row vector of length 'nVars' and return a
% scalar value. 
% -> 'nVars': Positive integer representing the number of variables
% in the problem. 
% -> 'lowerBound': Vector of lower bounds. 
% -> 'upperBound': Vector of upper bounds.
% -> *'constr' : Function handle. DE algorithm attempts to achieve
% constr(x) <=0.
% -> *'Options': Structure Array containing some of the folowing:
%   ---> 'populSize': Population size (default = 20);
%   ---> 'crossProb': Crossover probability (default = 0.6);
%   ---> 'diffWeight': Diferential Weight (default = 1);
%   ---> 'stopCriteria': Stop criteria based on diversity index
%   (default = 2e-5);
%   ---> 'maxGen':  Maximum number of generations (default = 100);
%   ---> 'penalityWeight': Penality weight (default = 100);
%   ---> 'penalityNumbOfGen': Number of generations 
%   adaptative penality weight scheme waits between actions (default = 3);
%   ---> 'penalityFactor': Factor by witch penality 'Weight' will
%   be multiplied or divided (default = 5);
%
% Outputs:
% -> 'xBest': Best point that DE algorithm located during its iterations.
% -> 'fBest': Fitness function evaluated at 'xBest'.
% -> 'nGen': Number of generations computed.
% -> 'X': Matrix cell containing all members of each
% interation. For example, X{1,2} is the second member of the first
% generation.
% -> 'fx': Matrix cell containing fitness function evaluated for
% each value of X
% -> 'gx': Matrix cell containing constraint evaluation for
% each value of X
%
% * optional parameters
%
% References:
% [1] Storn, R.; Price, K. (1997). "Differential evolution - a simple
% and efficient heuristic for global optimization over continuous
% spaces". Journal of Global Optimization 11:
% 341–359. doi:10.1023/A:1008202821328     
    


% Evaluate input arguments
% Check number of inputs
if (nargin > 6 || nargin < 4)
    error('requires at least 4 arguments and at most 5')
end
% Fill in unset optional values
switch nargin
  case 4
    constr = @(x) [];
    Options = struct;
  case 5
    Options = struct;
end

% Default Values
populSize = 7*nVars; % Population size
crossProb = 1; % Crossover probability (0.6,0.9)
diffWeight = 0.6; % Differential weight (0.4,1)
stopCriteria = 2e-5; % Stop criteria (based on diversity index)
maxGen = 10000; % Maximum number of generations



% fitness function
f = @(x) fitnessfcn(x);

% constraints
g = @(x) constraints(lowerBound, upperBound, constr, x);

nConstr =  length(g(zeros(nVars, 1))); % number of constraints
penalityWeight = 30; % penality weight
penalityNumbOfGen = [3*ones(nConstr,1)]; % adaptative penality # of generations
penalityFactor = [5*ones(nConstr,1)]; % adaptative penality factor

% Set optional values
names = fieldnames(Options);

for i = 1: length(names)
    switch char(names(i))
      case 'populSize'
        populSize = Options.populSize;
      case 'crossProb'
        crossProb = Options.crossProb;
      case 'diffWeight'
        diffWeight = Options.diffWeight;
      case 'stopCriteria'
        stopCriteria = Options.stopCriteria;
      case 'maxGen'
        maxGen = Options.maxGen;
      case 'penalityWeight'
        penalityWeight = Options.penalityWeight;
      case 'penalityNumbOfGen'
        penalityNumbOfGen = Options.genalityNumbOfGen;
      case 'penalityFactor'
        penalityFactor = Options.penalityFactor;
      otherwise
        error('unknow option');
    end
end



nGen = 1; % number of generations


% Initialize cells in withch fitness off all population will
% be store at each generation
fx = cell(maxGen,populSize);
gx = cell(maxGen,populSize);


% Initialize population (Random)
X = cell(maxGen,populSize);


maxF = -inf;
minF = inf;
sumG = zeros(nConstr,1);
nViolations = zeros(nConstr,1);


for i = 1:populSize
    % Initialize population (Random)
    for j = 1:nVars
        X{1,i}(j) = (upperBound(j) - lowerBound(j)).*rand + lowerBound(j);
    end
    
    % Get fitnessfcn and constraints violations
    fx{1,i} = f(X{1,i});
    gx{1,i} = g(X{1,i});
    
    % Gather information about maximum f
    if fx{1,i} > maxF
        maxF = fx{1,i};
    end
    
    % minimum f
    if fx{1,i} < minF
        minF = fx{1,i};
    end
    
    % sum of constraint violations
    for k = 1:nConstr
        if gx{1,i}(k)  > 0
            sumG(k) = sumG(k) + gx{1,i}(k);
            nViolations(k) = nViolations(k) + 1;
        end
    end
    
    
end

% Penality weight for each constraint
penalityWeight = [penalityWeight*(nViolations+1).*(maxF-minF)./(sumG+1)];



% Initialize mutant population (Zeros)

U = cell(populSize,1);
fu = cell(populSize,1);
gu = cell(populSize,1);

%% Start Algorithm
diversityIndex = stopCriteria + 1; % Diversity index used as stop criteria
genInConstr = zeros(nConstr,1); % Count number of generations
                                % xbest is respecting the i-th
                                % constraint. Reseting at each violation.
genOutConstr = zeros(nConstr,1); % Count number of generations
                                 % xbest is not respecting the i-th
                                 % constraint. Reseting at each violation.
while ~(nGen>maxGen || (diversityIndex < stopCriteria && sum(gBest) == 0))
    
    % initialize best individues of this generation
    fBest = inf;
    gBest = inf(nConstr,1);
    xBest = zeros(nVars,1);
    
    for i = 1:populSize
        
        % Ramdomly select r1,r2 e r3
        r=randperm(populSize,3);
        
        % Ramdomly select "di"
        di = randperm(nVars,1);
        
        % Crossover and mutation
        for j = 1:nVars
            if rand < crossProb || j == di
                U{i}(j) = X{nGen,r(1)}(j) + diffWeight.*(X{nGen,r(2)}(j) - X{nGen,r(3)}(j));
            else
                U{i}(j) = X{nGen,i}(j);
            end
        end
        
        % Calculate 'f' & 'g' 
        fu{i} = f(U{i}); gu{i} = g(U{i});
     end
        
     for i = 1:populSize
        % Selection
        if fu{i} + dot(penalityWeight,gu{i}) < fx{nGen,i} + dot(penalityWeight,gx{nGen,i})
            X{nGen+1,i} = U{i};
            fx{nGen+1,i} = fu{i};
            gx{nGen+1,i} = gu{i};
        else
            X{nGen+1,i} = X{nGen,i};
            fx{nGen+1,i} = fx{nGen,i};
            gx{nGen+1,i} = gx{nGen,i};
        end
        % Best individue
        if fBest + dot(penalityWeight,gBest) > fx{nGen,i}+ dot(penalityWeight,gx{nGen,i})
            xBest = X{nGen+1,i};
            fBest = fx{nGen+1};
            gBest = gx{nGen+1};
        end
    end
    

    % Reevaluate penality
    for i = 1:nConstr
        if gBest(i) == 0,
            genInConstr(i) = genInConstr(i) + 1;
            genOutConstr(i) = 0;
        else
            genOutConstr(i) = genOutConstr(i) + 1;
            genInConstr(i) = 0;
        end
        
        if genInConstr(i) == penalityNumbOfGen(i)-1,
            penalityWeight(i) = penalityWeight(i)/ penalityFactor(i);
        end
        
        if genOutConstr(i) == penalityNumbOfGen(i)-1,
            penalityWeight(i) = penalityWeight(i)* penalityFactor(i);
        end
    end
    
    % Stop Criteria
    mean = zeros(1,nVars);
    for i = 1:populSize
        mean =  X{nGen+1,i} + mean;
    end
    mean = mean/populSize;
    
    diversityIndex = 0;
    for i = 1:populSize
        diversityIndex =  norm(X{nGen+1,i} - mean)^2 + diversityIndex;
    end
    
    diversityIndex = diversityIndex/populSize;
    
    diversityIndex = diversityIndex/norm(upperBound - lowerBound);
    
    nGen = nGen+1;% Next generation
end
end

function g = constraints(lowerBound, upperBound, constr, x)
% CONSTRAINTS - Local Function: sumarize all
% constraints into one

i = 1;
    
for j = 1: length(lowerBound)
    g(i) = max(lowerBound(j) - x(j), 0);
    i = i+1;
end

for j = 1: length(upperBound)
    g(i) = max(x(j) - upperBound(j), 0);
    i = i+1;
end

C = constr(x);
for j = 1: length(C)
    g(i) = max(C(j), 0);
    i = i+1;
end

end
