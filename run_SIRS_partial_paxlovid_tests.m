clear all; close all;
par.pi_r = 1/(6/7); % five periods
par.pi = 2*par.pi_r ;
par.pi_s = 1/26; % twice a year
par.pi_lc = 0.04*0.08; % quite low
par.pi_lC1 = 0.8512362956 ;
par.pi_lC2 = 0.0497855435 ;
par.pi_lC3 = 0.0881757138 ;
par.piWt =  4/7;
par.piWu = 3/7;
par.share = 0.9;
par.test = [0.01 0.2 0.9];
[db,dbQ,dbA] = SIRSpaxtest(par,{'test = 0.01','test = 0.2','test = 0.9'});
dbA.Iw
dbA.LC
