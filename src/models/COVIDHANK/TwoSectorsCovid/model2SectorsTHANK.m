function [G1, C, impact, eu, SDX, zmat, NY, NX] = model2SectorsTHANK(param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For given parameter values, this function
% 1) puts the model of BDNPT in Gensys' canonical form
% 2) solves the RE system of equations using Chris Sims' Gensys
%
% The solution takes the form:  x(t) = G1 * x(t-1) + impact * e(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



pS          = 39;   % profits (S-type HHs)
% -------------------------------------------------------------------------
% INDEX for endogenous variables
% -------------------------------------------------------------------------
y           = 1;     % output
y1          = 2;     % output (sector 1)
y2          = 3;
k           = 4;     % capital
kbarS       = 5;     % physical capital (S-type HHs)
L1          = 6;     % hours (sector 1)
L2          = 7;     % hours (sector 2)
klr1        = 8;     % capital-labor ratio (sector 1)
klr2        = 9;     % capital-labor ratio (sector 2)
w1          = 10;     % real wage (sector 1)
w2          = 11;    % real wage (sector 2)
p           = 12;    % inflation
p1          = 13;    % inflation (sector 1)
p2          = 14;    % inflation (sector 2)
mc1         = 15;    % marginal cost (sector 1)
mc2         = 16;    % marginal cost (sector 2)
lambdap     = 17;    % price markup shock
lambdaw     = 18;    % wage markup shock
lambdaS     = 19;   % mu of consumption (S-type HHs)
lambdaH1    = 20;   % mu of consumption (H1-type HHs)
lambdaH2    = 21;   % mu of consumption (H2-type HHs)
lambdaSH1   = 22;   % average mu of consumption (S and H2-type HHs)
lambdaSH2   = 23;   % average mu of consumption (S and H2-type HHs)
c           = 24;   % consumption
cH1         = 25;   % consumption (H1-type HHs)
cH2         = 26;   % consumption (H2-type HHs)
cS          = 27;   % consumption (S-type HHs)
iS          = 28;   % investment (S-type HHs)
R           = 29;   % nominal interest rate
u           = 30;   % capital utilization
phi         = 31;   % "capital" multiplier
z           = 32;  % productivity shock
Rk          = 33;    % return on capital (rho)
a           = 34;  % relative technology
d           = 35;  % demand for sectoral goods
x           = 36;   % real GDP
b           = 37;  % intertemporal preference shock
t           = 38;   % taxes
g           = 39; % government spending shock
gw1         = 40;   % wage gap (sector 1)
gw2         = 41;   % wage gap (sector 1)
miu         = 42;  % investment shock
bR          = 43;   % real government debt
mp          = 44; % monetary policy shock
tS          = 45;   % taxes (S-type HHs)
tH1         = 46;   % taxes (H1-type HHs)
tH2         = 47;   % taxes (H2-type HHs)
tau         = 48;   % transfers
tauS        = 49; % transfer shock (S-type HHs)
tauH1       = 50;  % transfer shock (H1-type HHs)
tauH2       = 51;  % transfer shock (H2-type HHs)


% "expectation" variables
ep1         = 52;
ep2         = 53;
elambdaS    = 54;
elambdaH1   = 55;
elambdaH2   = 56;
ex          = 57;
ephi        = 58;
eRk         = 59;
eiS         = 60;
ew1         = 61;
ep          = 62;
ew2         = 63;

lastvar     = 63;

% exogenous processes



Nx          = lastvar + 12; % total number of endogenous variables


% -------------------------------------------------------------------------
% INDEX for Shocks
% -------------------------------------------------------------------------
as          = 1;
ds          = 2;
lambdaps    = 3;
bs          = 4;
zs          = 5;
tauH1s      = 6;
tauH2s      = 7;
mius        = 8;
lambdaws    = 9;
tauSs       = 10;
gs          = 11;
mps         = 12;

Nshocks     = 12;


% -------------------------------------------------------------------------
% INDEX for expectational errors
% -------------------------------------------------------------------------
p1ex       = 1;
p2ex       = 2;
lambdaSex  = 3;
lambdaH1ex = 4;
lambdaH2ex = 5;
xex        = 6;
phiex      = 7;
Rkex       = 8;
iSex       = 9;
w1ex       = 10;
pex        = 11;
w2ex       = 12;

Nex        = 12; % total number of expectational errors


% -------------------------------------------------------------------------
% parameter values and steady state
% -------------------------------------------------------------------------
SteadyState;

% put steady state values into a vector and save
ssvec = [gam disc 0 pss pSss exp(gam) L1ss L2ss Rkss w1ss w2ss ...
         klr1ss klr2ss kss F1 y1ss yss y2ss F2 xss iSss gss css ...
         tauss tauH1ss tauH2ss tauSss tss tH1ss tH2ss tSss cH10 ...
         cH20 cS0 lambdaH10 lambdaH20 lambdaS0 bR0 R0 mc1ss mc2ss];

% -------------------------------------------------------------------------
% System Matrices
% -------------------------------------------------------------------------
GAM0 = zeros(Nx, Nx);
GAM1 = zeros(Nx, Nx);
PSI  = zeros(Nx, Nshocks);
PPI  = zeros(Nx, Nex);
C    = zeros(Nx, 1);


% -------------------------------------------------------------------------
% equations
% -------------------------------------------------------------------------
GAM0(1,klr1) = -1;
GAM0(1,w1) = 1;
GAM0(1,Rk) = -1;

% -------------------------------------------------------------------------
GAM0(2,mc1) = -1;
GAM0(2,Rk) = alfa;
GAM0(2,w1) = 1-alfa;

% -------------------------------------------------------------------------
GAM0(3,klr2) = -1;
GAM0(3,w2) = 1;
GAM0(3,Rk) = -1;

% -------------------------------------------------------------------------
GAM0(4,mc2) = -1;
GAM0(4,Rk) = alfa;
GAM0(4,w2) = 1-alfa;
GAM0(4,a) = alfa-1;

% -------------------------------------------------------------------------
GAM0(5,y1) = -1;
GAM0(5,klr1) = (1+lambdapss)*alfa;
GAM0(5,L1) = (1+lambdapss);

% -------------------------------------------------------------------------
GAM0(6,y2) = -1;
GAM0(6,klr2) = (1+lambdapss)*alfa;
GAM0(6,L2) = (1+lambdapss);
GAM0(6,a) = -lambdapss;

% -------------------------------------------------------------------------
GAM0(7,y) = -1;
GAM0(7,y1) = v;
GAM0(7,y2) = 1-v;
GAM0(7,d) = (v-1)/(1-vars);

% -------------------------------------------------------------------------
GAM0(8,p) = -1;
GAM0(8,d) = (1-v)/(1-vars);
GAM0(8,p1) = v;
GAM0(8,p2) = 1-v;
GAM1(8,d) = (1-v)/(1-vars);

% -------------------------------------------------------------------------
GAM0(9,y1) = 1;
GAM0(9,y2) = -1;
GAM0(9,p1) = vars;
GAM0(9,p2) = -vars;
GAM0(9,d) = 1;
GAM1(9,y1) = 1;
GAM1(9,y2) = -1;
GAM1(9,d) = 1;

% -------------------------------------------------------------------------
GAM0(10,p1) = -1;
GAM0(10,ep1) = disc/(1+iotap*disc);
GAM0(10,mc1) = kappap;
GAM0(10,lambdap) = 1;
GAM1(10,p1) = -iotap/(1+iotap*disc);

% -------------------------------------------------------------------------
GAM0(11,p2) = -1;
GAM0(11,ep2) = disc/(1+iotap*disc);
GAM0(11,mc2) = kappap;
GAM0(11,lambdap) = 1;
GAM1(11,p2) = -iotap/(1+iotap*disc);

% -------------------------------------------------------------------------
GAM0(12,lambdaS) = -1;
GAM0(12,b) = 1;
GAM0(12,z) = -h/(exp(gam)-h);
GAM0(12,cS) = -exp(gam)/(exp(gam)-h);
GAM1(12,cS) = -h/(exp(gam)-h);

% -------------------------------------------------------------------------
GAM0(13,lambdaH1) = -1;
GAM0(13,b) = 1;
GAM0(13,z) = -h/(exp(gam)-h);
GAM0(13,cH1) = -exp(gam)/(exp(gam)-h);
GAM1(13,cH1) = -h/(exp(gam)-h);

% -------------------------------------------------------------------------
GAM0(14,lambdaH2) = -1;
GAM0(14,b) = 1;
GAM0(14,z) = -h/(exp(gam)-h);
GAM0(14,cH2) = -exp(gam)/(exp(gam)-h);
GAM1(14,cH2) = -h/(exp(gam)-h);

% -------------------------------------------------------------------------
GAM0(15,c) = -1;
GAM0(15,cS) = (1-theta)*cSss/css;
GAM0(15,cH1) = theta*fH1*cH1ss/css;
GAM0(15,cH2) = theta*(1-fH1)*cH2ss/css;

% -------------------------------------------------------------------------
GAM0(16,lambdaS) = -1;
GAM0(16,R) = 1;
GAM0(16,z) = -rhoz;
GAM0(16,ep) = -1;
GAM0(16,elambdaS) = disc*Rss*sss/exp(gam)/pss;
GAM0(16,elambdaH1) = disc*Rss*(1-sss)/exp(gam)/pss*fS1*lambdaH1ss/lambdaSss;
GAM0(16,elambdaH2) = disc*Rss*(1-sss)/exp(gam)/pss*(1-fS1)*lambdaH2ss/lambdaSss;
GAM0(16,ex) = disc*Rss*sss*sx/exp(gam)/pss*(1-fS1*lambdaH1ss/lambdaSss-(1-fS1)*lambdaH2ss/lambdaSss);

% -------------------------------------------------------------------------
GAM0(17,cH1) = -1;
GAM0(17,w1) = w1ss*H1ss/cH1ss;
GAM0(17,L1) = w1ss*H1ss/cH1ss;
GAM0(17,tH1) = -xss/cH1ss;
GAM0(17,tauH1) = xss/cH1ss;
GAM0(17,z) = -fS1*(1-sss)*Rss*bRss/(theta*fH1*exp(gam)*pss*cH1ss);
GAM0(17,p) = -fS1*(1-sss)*Rss*bRss/(theta*fH1*exp(gam)*pss*cH1ss);
GAM0(17,x) = -fS1*sss*sx*Rss*bRss/(theta*fH1*exp(gam)*pss*cH1ss);
GAM1(17,R) = -fS1*(1-sss)*Rss*bRss/(theta*fH1*exp(gam)*pss*cH1ss);
GAM1(17,bR) = -fS1*(1-sss)*Rss/(theta*fH1*exp(gam)*pss*cH1ss)*xss;

% -------------------------------------------------------------------------
GAM0(18,cH2) = -1;
GAM0(18,w2) = w2ss*H2ss/cH2ss;
GAM0(18,L2) = w2ss*H2ss/cH2ss;
GAM0(18,tH2) = -xss/cH2ss;
GAM0(18,tauH2) = xss/cH2ss;
GAM0(18,z) = -(1-fS1)*(1-sss)*Rss*bRss/(theta*(1-fH1)*exp(gam)*pss*cH1ss);
GAM0(18,p) = -(1-fS1)*(1-sss)*Rss*bRss/(theta*(1-fH1)*exp(gam)*pss*cH1ss);
GAM0(18,x) = -(1-fS1)*sss*sx*Rss*bRss/(theta*(1-fH1)*exp(gam)*pss*cH1ss);
GAM1(18,R) = -(1-fS1)*(1-sss)*Rss*bRss/(theta*(1-fH1)*exp(gam)*pss*cH1ss);
GAM1(18,bR) = -(1-fS1)*(1-sss)*Rss/(theta*(1-fH1)*exp(gam)*pss*cH2ss)*xss;

% -------------------------------------------------------------------------
GAM0(19,Rk) = 1;
GAM0(19,u) = -chi;

% -------------------------------------------------------------------------
GAM0(20,phi) = -1;
GAM0(20,ephi) = disc*exp(-gam)*(1-delta);;
GAM0(20,z) = -rhoz;
GAM0(20,elambdaS) = 1-disc*exp(-gam)*(1-delta);
GAM0(20,eRk) = 1-disc*exp(-gam)*(1-delta);

% -------------------------------------------------------------------------
GAM0(21,lambdaS) = -1;
GAM0(21,phi) = 1;
GAM0(21,miu) = 1;
GAM0(21,iS) = -exp(2*gam)*S*(1+disc);
GAM0(21,z) = -exp(2*gam)*S*(1-disc*rhoz);
GAM0(21,eiS) = exp(2*gam)*S*disc;
GAM1(21,iS) = -exp(2*gam)*S;

% -------------------------------------------------------------------------
GAM0(22,k) = -1;
GAM0(22,u) = 1;
GAM0(22,kbarS) = 1;
GAM0(22,z) = -1;

% -------------------------------------------------------------------------
GAM0(23,kbarS) = 1;
GAM0(23,z) = (1-delta)*exp(-gam);
GAM0(23,miu) = (1-delta)*exp(-gam)-1;
GAM0(23,iS) = (1-delta)*exp(-gam)-1;
GAM1(23,kbarS) = (1-delta)*exp(-gam);

% -------------------------------------------------------------------------
GAM0(24,w1) = -1;
GAM0(24,ew1) = disc/(1+disc);
GAM0(24,gw1) = -kappaw;
GAM0(24,p) = -(1+disc*iotaw)/(1+disc);
GAM0(24,ep) = disc/(1+disc);
GAM0(24,z) = -(1+disc*iotaw-rhoz*disc)/(1+disc);
GAM0(24,lambdaw) = 1;
GAM1(24,w1) = -1/(1+disc);
GAM1(24,p) = -iotaw/(1+disc);
GAM1(24,z) = -iotaw/(1+disc);

% -------------------------------------------------------------------------
GAM0(25,gw1) = 1;
GAM0(25,w1) = -1;
GAM0(25,L1) = niu;
GAM0(25,b) = 1;
GAM0(25,lambdaSH1) = -1;

% -------------------------------------------------------------------------
GAM0(26,lambdaSH1) = -1;
GAM0(26,lambdaS) = (1-theta)*fS1*lambdaSss/[(1-theta)*fS1+theta*fH1]/lambdaSH1ss;
GAM0(26,lambdaH1) = theta*fH1*lambdaH1ss/[(1-theta)*fS1+theta*fH1]/lambdaSH1ss;

% -------------------------------------------------------------------------
GAM0(27,w2) = -1;
GAM0(27,ew2) = disc/(1+disc);
GAM0(27,gw2) = -kappaw;
GAM0(27,p) = -(1+disc*iotaw)/(1+disc);
GAM0(27,ep) = disc/(1+disc);
GAM0(27,z) = -(1+disc*iotaw-rhoz*disc)/(1+disc);
GAM0(27,lambdaw) = 1;
GAM1(27,w2) = -1/(1+disc);
GAM1(27,p) = -iotaw/(1+disc);
GAM1(27,z) = -iotaw/(1+disc);

% -------------------------------------------------------------------------
GAM0(28,gw2) = 1;
GAM0(28,w2) = -1;
GAM0(28,L2) = niu;
GAM0(28,b) = 1;
GAM0(28,lambdaSH2) = -1;

% -------------------------------------------------------------------------
GAM0(29,lambdaSH2) = -1;
GAM0(29,lambdaS) = (1-theta)*(1-fS1)*lambdaSss/[(1-theta)*(1-fS1)+theta*(1-fH1)]/lambdaSH2ss;
GAM0(29,lambdaH2) = theta*(1-fS1)*lambdaH2ss/[(1-theta)*(1-fS1)+theta*(1-fH1)]/lambdaSH2ss;

% -------------------------------------------------------------------------
GAM0(30,R) = -1;
GAM0(30,p) = (1-rhoR)*fp;
GAM0(30,x) = (1-rhoR)*(fx+fdx);
GAM0(30,mp) = 1;
GAM1(30,R) = -rhoR;
GAM1(30,x) = (1-rhoR)*fdx;

% -------------------------------------------------------------------------
GAM0(31,t) = -1;
GAM0(31,tS) = 1-theta;
GAM0(31,tH1) = theta*fH1;
GAM0(31,tH2) = theta*(1-fH1);

% -------------------------------------------------------------------------
GAM0(32,tau) = -1;
GAM0(32,tauS) = 1-theta;
GAM0(32,tauH1) = theta*fH1;
GAM0(32,tauH2) = theta*(1-fH1);

% -------------------------------------------------------------------------
% GAM0(33,bR) = xss/bRss;
% GAM0(33,p) = Rss/exp(gam)/pss;
% GAM0(33,z) = Rss/exp(gam)/pss;
% GAM0(33,g) = -(1-Rss/exp(gam)/pss)*xss/(gss-tss+tauss);
% GAM0(33,t) = (1-Rss/exp(gam)/pss)*xss/(gss-tss+tauss);
% GAM0(33,tau) = -(1-Rss/exp(gam)/pss)*xss/(gss-tss+tauss);
% GAM1(33,R) = Rss/exp(gam)/pss;
% GAM1(33,bR) = Rss/exp(gam)/pss*xss/bRss;

GAM0(33,bR) = xss*(gss-tss+tauss);
GAM0(33,p) = Rss/exp(gam)/pss*bRss*(gss-tss+tauss);
GAM0(33,z) = Rss/exp(gam)/pss*bRss*(gss-tss+tauss);
GAM0(33,g) = -(1-Rss/exp(gam)/pss)*xss*bRss;
GAM0(33,t) = (1-Rss/exp(gam)/pss)*xss*bRss;
GAM0(33,tau) = -(1-Rss/exp(gam)/pss)*xss*bRss;
GAM1(33,R) = Rss/exp(gam)/pss*bRss*(gss-tss+tauss);
GAM1(33,bR) = Rss/exp(gam)/pss*xss*(gss-tss+tauss);

% -------------------------------------------------------------------------
GAM0(34,t) = -xss/tss;
GAM0(34,x) = 1+(1-rhot)*ziX-(1-rhot)*ziG;
GAM0(34,g) = (1-rhot)*ziG*xss/gss;
GAM1(34,t) = -rhot*xss/tss;
GAM1(34,bR) = -(1-rhot)*ziB*xss/bRss;
GAM1(34,x) = rhot+(1-rhot)*ziB;

% -------------------------------------------------------------------------
GAM0(35,tH1) = theta*fH1;
GAM0(35,t) = -psiH1;

% -------------------------------------------------------------------------
GAM0(36,tH2) = theta*(1-fH1);
GAM0(36,t) = -psiH2;

% -------------------------------------------------------------------------
GAM0(37,pS) = -1;
GAM0(37,y) = yss/(1-theta)/pSss;
GAM0(37,mc1) = -mc1ss*klr1ss*L1ss/(1-theta)/pSss;
GAM0(37,klr1) = -mc1ss*klr1ss*L1ss/(1-theta)/pSss;
GAM0(37,L1) = -mc1ss*klr1ss*L1ss/(1-theta)/pSss;
GAM0(37,mc2) = -A2^(1-alfa)*mc2ss*klr2ss*L2ss/(1-theta)/pSss;
GAM0(37,klr2) = -A2^(1-alfa)*mc2ss*klr2ss*L2ss/(1-theta)/pSss;
GAM0(37,L2) = -A2^(1-alfa)*mc2ss*klr2ss*L2ss/(1-theta)/pSss;

% -------------------------------------------------------------------------
GAM0(38,klr1) = klr1ss*L1ss/kss;
GAM0(38,L1) = klr1ss*L1ss/kss;
GAM0(38,klr2) = klr2ss*L2ss/kss;
GAM0(38,L2) = klr2ss*L2ss/kss;
GAM0(38,k) = -1;

% -------------------------------------------------------------------------
GAM0(39,c) = css/yss;
GAM0(39,iS) = (1-theta)*iSss/yss;
GAM0(39,g) = xss/yss;
GAM0(39,u) = Rkss*kss/yss;
GAM0(39,y) = -1;

% -------------------------------------------------------------------------
GAM0(40,x) = -1;
GAM0(40,y) = 1;
GAM0(40,u) = -Rkss*kss/yss;

% -------------------------------------------------------------------------
GAM0(41,a) = 1;         GAM1(41,a) = rhoa;              PSI(41,as) = 1;
GAM0(42,d) = 1;         GAM1(42,d) = rhod;              PSI(42,ds) = 1;
GAM0(43,lambdap) = 1;   GAM1(43,lambdap) = rholambdap;  PSI(43,lambdaps) = 1;
GAM0(44,b) = 1;         GAM1(44,d) = rhob;              PSI(44,bs) = 1;
GAM0(45,z) = 1;         GAM1(45,z) = rhoz;              PSI(45,zs) = 1;
GAM0(46,tauH1) = 1;     GAM1(46,tauH1) = rhotauH1;      PSI(46,tauH1s) = 1;
GAM0(47,tauH2) = 1;     GAM1(47,tauH2) = rhotauH2;      PSI(47,tauH2s) = 1;
GAM0(48,miu) = 1;       GAM1(48,miu) = rhomiu;          PSI(48,mius) = 1;
GAM0(49,lambdaw) = 1;   GAM1(49,lambdaw) = rholambdaw;  PSI(49,lambdaws) = 1;
GAM0(50,tauS) = 1;      GAM1(50,tauS) = rhotauS;        PSI(50,tauSs) = 1;
GAM0(51,g) = 1;         GAM1(51,g) = rhog;              PSI(51,gs) = 1;
GAM0(52,mp) = 1;        GAM1(52,mp) = rhomp;            PSI(52,mp) = 1;

% -------------------------------------------------------------------------
GAM0(53,p1) = 1; GAM1(53,ep1) = 1; PPI(53,p1ex) = 1;
GAM0(54,p2) = 1; GAM1(54,ep2) = 1; PPI(54,p2ex) = 1;
GAM0(55,lambdaS) = 1; GAM1(55,elambdaS) = 1; PPI(55,lambdaSex) = 1;
GAM0(56,lambdaH1) = 1; GAM1(56,elambdaH1) = 1; PPI(56,lambdaH1ex) = 1;
GAM0(57,lambdaH2) = 1; GAM1(57,elambdaH2) = 1; PPI(57,lambdaH2ex) = 1;
GAM0(58,x) = 1; GAM1(58,ex) = 1; PPI(58,xex) = 1;
GAM0(59,phi) = 1; GAM1(59,ephi) = 1; PPI(59,phiex) = 1;
GAM0(60,Rk) = 1; GAM1(60,eRk) = 1; PPI(60,Rkex) = 1;
GAM0(61,iS) = 1; GAM1(61,eiS) = 1; PPI(61,iSex) = 1;
GAM0(62,w1) = 1; GAM1(62,ew1) = 1; PPI(62,w1ex) = 1;
GAM0(63,p) = 1; GAM1(63,ep) = 1; PPI(63,pex) = 1;
GAM0(64,w2) = 1; GAM1(64,ew2) = 1; PPI(64,w2ex) = 1;

save('eqcond.mat', 'ssvec','GAM0','GAM1', 'C', 'PSI', 'PI');
% -------------------------------------------------------------------------
% Solution of the RE system of equations using Chris Sims' Gensys
% -------------------------------------------------------------------------
[G1, C, impact, fmat, fwt, ywt, gev, eu] = GENSYS(GAM0, GAM1, C, PSI, PPI);
