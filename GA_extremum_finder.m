function [ output_args ] = GA_extremum_finder( fun, n, bitstringlength, pc, pm, boundaries)

%   Takes a function at an input and attemps to 
%   determine either the minimum or the maxiumum.

B = bitstringlength;

%generate initial guess
pop = (boundaries(end) - boundaries(1))*rand(12,1);

%eval fitness
fitness = zeros(size(pop));
for i = 1:length(pop)
    fitness(i) = fun(pop(i));
end


init = dec2bin(pop,B);
output_args = fitness;

end

