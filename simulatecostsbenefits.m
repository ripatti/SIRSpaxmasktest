function [costs,benefits] = simulatecostsbenefits(m,par,qrange,qGDP,N,paxlovidcost,paxlovidshare,vaccinecost,vaccineshare,maskcost, maskshare, testcost, testshare)
% simulates SIRS and economic model for policy mix given by policy
% structure
% Inputs:
%  m model structure for the economics model
%  qrange simulation range
%  qGDP quarterly gdp volume in billions 
%  N population in millions
%  policy is a structure of policy mixes
% Ouputs:
%  cost  costs in bill USD
%  benefits benefits in bill US
[db,dbQ,dbA] = SIRSpaxtest(par,{'baseline','policy'});
g = simulateLaborSupply(m,qrange,dbQ.Iw+dbQ.LC);
gdp = g.y;
gdpA = 100*convert(gdp,'A','method=',@sum);
qCosts = qGDP*g.y;
aCosts = convert(qCosts, 'A', 'method=',@sum);
% costs of paxlovid
uspop =  N*1e6; % 18 and over, 2023, in https://data.census.gov/table/ACSDP1Y2023.DP05?q=United%20States&table=DP05&g=010XX00US&lastDisplayedRow=29&vintage=2017&layer=state&cid=DP05_0001E
%paxlovidcosts = db.NewCases*uspop*530/1000000000; % in billions of $
paxlovidcosts = db.NewCases*uspop*paxlovidcost*paxlovidshare/1000000000; % in billions of $
paxCostsA = convert(paxlovidcosts,'A', 'method=',@sum);
% I think here should be quarterly costs
maskCostsA = convert(uspop*2*(365/4)*maskcost*maskshare*Series(qrange,1)/1e9,'A', 'method=',@sum);
vaccinecostsA = convert(uspop*1*vaccinecost*vaccineshare*Series(qrange,1)/1e9,'A', 'method=',@sum);
testCostsA = convert(uspop*1*(52/4)*testcost*testshare*Series(qrange,1)/1e9,'A', 'method=',@sum);
benefits = aCosts;
benefits = -benefits{:,1} + benefits{:,2};
costs = paxCostsA + maskCostsA + vaccinecostsA + testCostsA;
costs = costs{:,2};
benefits(1)
costs(1)



