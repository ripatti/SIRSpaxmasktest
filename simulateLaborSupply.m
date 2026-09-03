function db = simulateLaborSupply(m,rng,sh)
% Simulates labor supply shock
% Input:    m       model
%           rng     simulation range
%           sh      shock time series
%           d
%
% (c) Antti Ripatti, 2024- 
p1 = Plan.forModel(m, rng);
p1 = exogenize(p1, rng, "myShock");
p1 = endogenize(p1, rng, "e_myShock");
d = zerodb(m,rng);
d1 = d;
d1.myShock = log(1-sh);
d1.myShock(rng(1)-1) = log(1-sh(rng(1))) ;
db = simulate( ...
  m, d1, rng ...
  , "plan", p1 ...
  ); 
