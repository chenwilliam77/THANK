clear all;
addpath ../../../util;
addpath ../../../util/gensys;
lambdapss   = 0.15;                 % net markup
alfa        = 1-.60/(1+lambdapss);  % to match CE/(GDP-PI)~0.6 (or 0.58?)
v           = .7;
vars        = 1.5;
disc        = .99;                  % implies annual 6% return on capital
xip         = .85;
iotap       = .21;
h           = .5;
gam         = log(1.02)/4;          % 2% annual growth rate
theta       = 1/4;
fH1         = 1/4;
pss         = 1.02^(1/4);           % 2% annualized inflation
sss         = .995;
fS1         = 3/4;
sx          = 0;
chi         = 5;
delta       = .025;                 % to match (I+Cdurables)/GDP=0.025
S           = 2.5;
iotaw       = .15;
xiw         = .75;
niu         = 4;
lambdawss   = .12;
fp          = 2;
fx          = .1;
fdx         = .25;
g_x         = .18;
t_x         = .2;
tau_x       = 0;
ziX         = 0;
ziG         = 0;
ziB         = 1;
psiH1       = 0;
psiH2       = 0;
A2          = .7;
H1ss        = 1;
H2ss        = 1;
rhoR        = .8;
rhot        = .8;

rhoa        = .15;
rhod        = .15;
rholambdap  = .95;
rhob        = .7;
rhoz        = .25;
rhotauH1    = .98;
rhotauH2    = .98;
rhomiu      = .7;
rholambdaw  = .98;
rhotauS     = .98;
rhog        = .98;
rhomp       = .8;

param=[lambdapss alfa v vars disc xip iotap h gam theta fH1 pss sss fS1 sx...
        chi delta S iotaw xiw niu lambdawss fp fx fdx g_x t_x tau_x ziX ziG...
        ziB  psiH1 psiH2 A2 H1ss H2ss rhoR rhot rhoa rhod rholambdap...
        rhob rhoz rhotauH1 rhotauH2 rhomiu rholambdaw rhotauS rhog rhomp];

[G1, C, impact, eu, SDX, zmat, NY, NX] = model2SectorsTHANK(param);
