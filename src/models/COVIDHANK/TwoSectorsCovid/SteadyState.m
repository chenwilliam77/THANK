% -------------------------------------------------------------------------
% INDEX for the parameters
% -------------------------------------------------------------------------
lambdapss   = param(1);
alfa        = param(2);
v           = param(3);
vars        = param(4);
disc        = param(5);
xip         = param(6);
iotap       = param(7);
h           = param(8);
gam         = param(9);
theta       = param(10);
fH1         = param(11);
pss         = param(12);
sss         = param(13);
fS1         = param(14);
sx          = param(15);
chi         = param(16);
delta       = param(17);
S           = param(18);
iotaw       = param(19);
xiw         = param(20);
niu         = param(21);
lambdawss   = param(22);
fp          = param(23);
fx          = param(24);
fdx         = param(25);
g_x         = param(26);
t_x         = param(27);
tau_x       = param(28);
ziX         = param(29);
ziG         = param(30);
ziB         = param(31);
psiH1       = param(32);
psiH2       = param(33);
A2          = param(34);
H1ss        = param(35);
H2ss        = param(36);
rhoR        = param(37);
rhot        = param(38);

rhoa        = param(39);
rhod        = param(40);
rholambdap  = param(41);
rhob        = param(42);
rhoz        = param(43);
rhotauH1    = param(44);
rhotauH2    = param(45);
rhomiu      = param(46);
rholambdaw  = param(47);
rhotauS     = param(48);
rhog        = param(49);
rhomp       = param(50);

kappap      = (1-xip*disc)*(1-xip)/xip/(1+iotap*disc);
kappaw      = (1-xiw*disc)*(1-xiw)/xiw/(1+disc)/(1+niu*(1+1/lambdawss));

tauH1ss=0;
tauH2ss=0;
tauSss=0;


% -------------------------------------------------------------------------
% steady state computation
% -------------------------------------------------------------------------
L1ss=[(1-theta)*fS1+theta*fH1]*H1ss;
L2ss=[(1-theta)*(1-fS1)+theta*(1-fH1)]*H2ss;

Rkss=exp(gam)/disc-(1-delta);
w1ss=[alfa^alfa*(1-alfa)^(1-alfa)/(1+lambdapss)/(Rkss^alfa)]^(1/(1-alfa));
w2ss=A2*w1ss;

klr1ss=w1ss/Rkss*alfa/(1-alfa);
klr2ss=A2*klr1ss;
kss=klr1ss*L1ss+klr2ss*L2ss;
F1=lambdapss/(1+lambdapss)*L1ss*klr1ss^alfa;
y1ss=L1ss*klr1ss^alfa-F1;
yss=y1ss/v;
y2ss=(1-v)*yss;
F2=lambdapss*y2ss/A2;
xss=yss;
iSss=[1-(1-delta)*exp(-gam)]*exp(gam)*kss/(1-theta);
gss=g_x*xss;
css=yss-gss-(1-theta)*iSss;

mc1ss=Rkss^alfa*w1ss^(1-alfa)/[alfa^alfa*(1-alfa)^(1-alfa)];
mc2ss=Rkss^alfa*(w2ss/A2)^(1-alfa)/[alfa^alfa*(1-alfa)^(1-alfa)];
pSss=[yss-mc1ss*L1ss*klr1ss^alfa-mc2ss*L2ss*klr2ss^alfa*A2^(1-alfa)]/(1-theta);

tauss=tau_x*xss;
tauH1ss=tauss*psiH1/theta/fH1;
tauH2ss=tauss*psiH2/theta/(1-fH1);
tauSss=(tauss-theta*(fH1*tauH1ss+(1-fH1)*tauH2ss))/(1-theta);

tss=t_x*xss;
tH1ss=tss*psiH1/theta/fH1;
tH2ss=tss*psiH2/theta/(1-fH1);
tSss=(tss-theta*(fH1*tH1ss+(1-fH1)*tH2ss))/(1-theta);

global w1ss w2ss H1ss H2ss tH1ss tH2ss tauH1ss tauH2ss css gss tss tauss
global sss theta fH1 fS1 gam pss h disc
cH10=w1ss*H1ss-tH1ss+tauH1ss;
cH20=w2ss*H2ss-tH2ss+tauH2ss;
cS0=(css-cH10*theta*fH1-cH20*theta*(1-fH1))/(1-theta);
lambdaH10=exp(gam)/(exp(gam)-h)/cH10;
lambdaH20=exp(gam)/(exp(gam)-h)/cH20;
lambdaS0=exp(gam)/(exp(gam)-h)/cS0;
bR0=(gss-tss+tauss)/(1-(1+.025/4)/exp(gam)/pss);
R0=(1+.025/4);

% ,'Display','iter'
options = optimoptions('fsolve','MaxFunctionEvaluations',8000,'MaxIterations',4000,'Display','final-detailed');
fun=@(x) distSS(x);
ss=fsolve(fun,[cH10 cH20 cS0 lambdaH10 lambdaH20 lambdaS0 bR0 R0],options)
cH1ss=ss(1);
cH2ss=ss(2);
cSss=ss(3);
lambdaH1ss=ss(4);
lambdaH2ss=ss(5);
lambdaSss=ss(6);
bRss=ss(7);
Rss=ss(8);

lambdaSH1ss=(1-theta)*fS1/[(1-theta)*fS1+theta*fH1]*lambdaSss+theta*fH1/[(1-theta)*fS1+theta*fH1]*lambdaH1ss;
lambdaSH2ss=(1-theta)*(1-fS1)/[(1-theta)*(1-fS1)+theta*(1-fH1)]*lambdaSss+theta*(1-fH1)/[(1-theta)*(1-fS1)+theta*(1-fH1)]*lambdaH2ss;


function r=distSS(x)
global w1ss w2ss H1ss H2ss tH1ss tH2ss tauH1ss tauH2ss css gss tss tauss
global sss theta fH1 fS1 gam pss h disc

cH1ss=x(1);
cH2ss=x(2);
cSss=x(3);
lambdaH1ss=x(4);
lambdaH2ss=x(5);
lambdaSss=x(6);
bRss=x(7);
Rss=x(8);

r(1)=w1ss*H1ss-tH1ss+tauH1ss+Rss*bRss*fS1*(1-sss)/theta/fH1/exp(gam)/pss-cH1ss;
r(2)=w2ss*H2ss-tH2ss+tauH2ss+Rss*bRss*(1-fS1)*(1-sss)/theta/(1-fH1)/exp(gam)/pss-cH2ss;
r(3)=css/(1-theta)-theta*(fH1*cH1ss+(1-fH1)*cH2ss)/(1-theta)-cSss;

r(4)=exp(gam)/(exp(gam)-h)/cH1ss-lambdaH1ss;
r(5)=exp(gam)/(exp(gam)-h)/cH2ss-lambdaH2ss;
r(6)=exp(gam)/(exp(gam)-h)/cSss-lambdaSss;

r(7)=(gss-tss+tauss)/(1-Rss/exp(gam)/pss)-bRss;

r(8)=exp(gam)*pss/disc*lambdaSss/[sss*lambdaSss+(1-sss)*(fS1*lambdaH1ss+(1-fS1)*lambdaH2ss)]-Rss;

end
