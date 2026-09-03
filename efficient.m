% efficient.m
% Simulating efficient frontier
% The policy tools: 
% 1) Paxlovid
% 2) Vaccines
% 3) Masks
% 4) Rapid tests
%
%% Clear Workspace
%
% Clear workspace, close all graphics figures, clear command window, and
% check the IRIS version.
clear all; close all; clc;
%% Control what to plot and what are the shocks
beta = 1/(1.03);
irfLength = 40;
plotList = { '"Hours" lab', '"Capital stock" k', '"rental rate" rk', '"Output" y' , '"Consumption" c', '"Investments" inve', ...
    '"Inflation (annualised)", pinfann', '"Nominal wage" w',  ...
    '"Interest rate (annualised)"  rann',  '"Labor supply shock" myShock'};

plotSlideList = {  '"Hours" lab',  '"Output" y' , '"consumption" c', '"Nomincal wages" w', '"Inflation (annualised)", pinfann', '"Interest rate (annualised)" rann'};
%% Read model
fprintf(1,'Preparing model...');
flags.isonlywedge = true(); % false = do not consider myShock in full fledge manner
m = Model.fromFile('SWd.model');
m = refresh(m);
m = solve(m); % IRIS8
flags.isonlywedge = false();
m1 = Model.fromFile('SWls.model',assign=flags);
m1 = refresh(m1);
m1 = solve(m1); % IRIS8
disp('Ready');
%% Cost parameters
qGDP = 29349.9/4; % GDO 2024q3 in USD billion
USpop =  262.266460;; % US adult popultaion million people
Cpaxlovidtreatment = 280; % costs dollars per dose
Cmasks = 1; % $ 1 per person per day
Cvaccine = 50; % $ 50 per person, entire population
Ctests = 0.5; % dollars per person, weekly testing test -share of entire population
%% Baseline parameters
% Assing parameter values
par.pi_ru = 1/(7.33/7);
par.pi_rt = 1/(3.73/7);
par.share = 0;
par.piu = 2*par.pi_ru; % =R0*pi_r
par.pit = 1.27*par.pi_rt;
% Let (1-pi_s) = P(staying immune), then average time of staying immune is
%             1/[1 - (1-pi_s)]
% Assume that the length of immunity is 4 months (16 weeks). Then
% 1/[1 - (1-pi_s)] = 16 and pi_s = 1/16.
par.pi_s = 1/26; % 1/(length of immunity in weeks), 1/26 half year frequency
par.pi_lcu = 0.04*0.08; % share of serious LC cases from new cases: incidence x 
par.pi_lct = 0.74*0.04*0.08;
%  pi_lclen = (1-0.2)^(1/52); % 20 percent still ill after one year
% pi_lclen = 0.662; % estimate from persistence of men tiredness from NL data https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)01214-4/fulltext#fig2
par.pi_lC1 = 0.8512362956 ;
par.pi_lC2 = 0.0497855435 ;
par.pi_lC3 = 0.0881757138 ;
par.piWt = 4/7;
par.piWu = 3/7;
par.test = 0.0;
[db,dbQ,dbA] = SIRSpaxtest(par,{'baseline'});
qrange = get(dbQ.S, 'range');
qrange = qrange(1):qq(2030,4);
plotRng = qrange(1):qq(2030,4);
%% Just loop over alternatives
Arange = yy(get(qrange(1),'year')):yy(get(qrange(end),'year'));
sharestep = 0.1;
i = 1; j=1; k=1;
for share = 0:sharestep:1
  par.share = [0 share];
  j = 1;
  for vaccshare = 0:sharestep:1
    par.piu = [2*par.pi_ru (1-vaccshare*0.4)*2*par.pi_ru]; 
    par.pit = [1.27*par.pi_rt (1-vaccshare*0.4)*1.27*par.pi_rt]; 
    k = 1;
    for maskshare = 0:sharestep:1
      par.piu = [par.piu(1) (1-maskshare*0.9)*par.piu(2)];
      par.pit = [par.pit(1) (1-maskshare*0.9)*par.pit(2)];
      l = 1;
      for testshare = 0:1:1
        par.test = testshare;
        [costsTS,benefitsTS] = simulatecostsbenefits(m,par,qrange,qGDP,USpop,Cpaxlovidtreatment,share,Cvaccine,vaccshare,Cmasks, maskshare, Ctests, testshare);
        costs(i,j,k) = sum(beta.^(0:numel(Arange)-1)' .* costsTS(Arange));
        benefits(i,j,k) = sum(beta.^(0:numel(Arange)-1)' .* benefitsTS(Arange));
        policy{i,j,k} = sprintf('(%3.2f, %3.2f, %3.2f, %3.2f)',share,vaccshare,maskshare,testshare);
        fprintf('.');
        l = l + 1;
      end
      k = k+1;
    end
    j = j+1;
    fprintf('\n');
  end
  i = i+1;
end
%% Creating figures
% efficient frontier
origCosts = costs;
origBenefits = benefits;
origPolicy = policy;
benefits = benefits(:);
costs = costs(:);
policy = policy(:);   % string array or cell array of char

isEfficient = true(size(benefits));

for i = 1:numel(benefits)
    isEfficient(i) = ~any( ...
        benefits >= benefits(i) & ...
        costs <= costs(i) & ...
        (benefits > benefits(i) | costs < costs(i)) ...
    );
end

% Efficient points and corresponding policy labels
benefitsEff = benefits(isEfficient);
costsEff = costs(isEfficient);
policyEff = policy(isEfficient);

% Sort frontier by benefit
[benefitsEff, order] = sort(benefitsEff);
costsEff = costsEff(order);
policyEff = policyEff(order);

figure;
scatter(benefits, costs, 40, 'filled');
hold on;
plot(benefitsEff, costsEff, '-o', 'LineWidth', 2);
scatter(benefitsEff, costsEff, 80, 'filled');
% Label only frontier points
text(benefitsEff, costsEff, policyEff, ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'left');
xlabel('Benefits');
ylabel('Costs');
title('Efficient frontier: costs vs benefits: (Paxlovid, vaccination, mask, test) shares');
grid on;
legend('All points', 'Efficient frontier', 'Efficient points', ...
    'Location', 'best');
hold off;
orient landscape
print -dpdf -fillpage figs\costsbenefits.pdf
%% Save simulation data
save('efficientfrontier.mat','origCosts','origBenefits','origPolicy','costs','benefits','policy');
