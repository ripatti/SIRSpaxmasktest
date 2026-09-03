function [db,dbQ,dbA] =  SIRSpaxtest(par,legendtxt)
% run SIR model with parameter values defined in input par
% running the SIR model with Iris toobox
%
% (c) Antti Ripatti, 2026-
%
% parameter values, data ranges
vacc = 1; % effective vaccination rate
range = ww(2022,40):ww(2052,52);
range = ww(2022,1):ww(2052,52);
initRange = ww(2021,52-36):ww(2052,52); % range for initial values
plotRange = range(1):ww(2024,52);
% read model
m = Model.fromFile('SIRSpaxtest.model','linear', false,'growth',false);
% assign supplied parameter values to the model
m = alter(m,max(struct2array((structfun(@size,par,'UniformOutput',false))))); % check the max dimension of parameter values
m = assign(m,par);
%% input database of initial values
%   S, lam, D, Du, Dt, I,  Iu, It, Iwt, Iwu, Iw, R, LC, NewCases
idb = struct();
idb.Iu = Series(initRange,0.01);
idb.It = Series(initRange,0.01);
idb.I = idb.Iu+idb.It;
idb.lam = Series(initRange,0);
idb.S = vacc-idb.I;
idb.R = 1 - idb.S - idb.I;
idb.NewCases = 0*idb.I;
idb.Iwt = par.piWt*idb.It;
idb.Iwu = par.piWu*idb.Iu;
idb.Iw = idb.Iwt + idb.Iwu;
idb.D = Series(initRange,0);
idb.Du = Series(initRange,0);
idb.Dt = Series(initRange,0);

% Exogenous arrival of new variants based on US wastewater data, https://covid.cdc.gov/covid-data-tracker/#wastewater-surveillance
idb.e_variant = Series(initRange,0);
% idb.e_variant(ww(2022,15)) = 0.2; % Apr 11, 2022
% idb.e_variant(ww(2022,46)) = 0.1; % Nov 19 11, 2022
% %idb.e_variant(ww(2023,29)) = 0.01; % Apr 11, 2023
% idb.e_variant(ww(2023,45)) = 0.3; % Nov 11, 2023
% idb.e_variant(ww(2024,32):ww(2024,40)) = 0.02; % Late summer peak 2024
len = length(initRange);
trendvals = 0.5:(0.5/26):1;
lentv = length(trendvals);
vals = [trendvals ones(1,len-lentv)];
idb.LC = 0.05*(1/12)*Series(initRange,vals); % see XLS file
%idb.LC = 0.06*0.15*Series(range,1); % see XLS file
idb.e_vacc = Series(range,0);
% p.pi =  0.1:0.05:0.2;
% m = assign(m, p);
%[db, ExitFlag, AddF, Delta]
solverOpts = {"iris-qnsd", "maxIterations", 100000, "functionTolerance", 1e-4,"display","notify"}; 
db = simulate(m, idb,range(1:end),'method','period', 'startIterationsfrom','data', 'solver', solverOpts); %,'solve','@lsqnonlin'); %,'nonlinear=',len);
%db = simulate(m, idb,1:len,'method','stacked','solve','@lsqnonlin'); %,'nonlinear=',len);
%db.Re = m.pi/m.pi_r*db.S;
db.LCpop = 5500000.*db.LC;
%db.NewCases(range(1)) = 0;
%db.R(1:len,2) = (1-1/R0)*ones(len,1); 
%db.S(1:len,2) = (1/R0)*ones(len,1); 
%db.Re(1:len,2) = ones(len,1); 
%{
figure();
plot([ db.I db.LC])
legend('Infected','long COVID');
figure;
dbplot(db,ww(2022,1):ww(2030,52),{'"Susceptible" S','"Infected (treated)" It','"Infected (untreated)" Iu','"Recovered" R','"Serious long COVID" LC'});
grfun.bottomlegend(legendtxt);
%legend('Location', 'southoutside', 'Orientation', 'horizontal');
orient landscape
print -dpdf -fillpage SIRSpaxlovid.pdf
%}
%% temporal aggregation
dbQ = struct();
dbA = struct();
fn = fieldnames(db);
for i=1:length(fn)
  if length(db.(fn{i})) > 1
    dbQ.(fn{i}) = convert(db.(fn{i}),'Q');
    dbA.(fn{i}) = convert(db.(fn{i}),'A');
  end
end
% shorter sick leaves asymptomatic infections
save('dbQ.mat','dbQ');
%% Steady state I and LC and storing the results
% endemicSS(m.pi_r,R0,pi_s,pi_lC,rho)
%[ss.I,ss.LC] = endemicSS(m.pi_r,m.pi./m.pi_r,m.pi_s,m.pi_lc,m.pi_lC1, m.pi_lC2, m.pi_lC3, vacc);
save('dbQ.mat','dbQ', 'dbA','db','m');

