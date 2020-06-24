%function ss=SteadyStateCovid(par);

lambdap=0.15;               % net markup
alpha=1-.60/(1+lambdap);    % to match CE/(GDP-PI)~0.6 (or 0.58?)
delta=.025;                 % to match (I+Cdurables)/GDP=0.25
piss=1.02^(1/4);            % 2% annualized inflation
gam=log(1.02)/4;            % 2% annual growth rate
disc=.99;                   % implies annual 6% return on capital

h=.5;                       % habit formation (really important for MU)
g_x=.18;                    % to match G/(GDP - net exports) = 0.18

s=.995;
theta=1/4;
fs1=3/4;
fh1=1/4;
A2=.7;
v=.7;
t_x=.2;
tau_x=0;
%psih1=theta*fh1/1.1;
%psih2=theta*(1-fh1)/1.1;
psih1 = 0;
psih2 = 0;

H1=1;
H2=1;


L1=((1-theta)*fs1+theta*fh1)*H1;
L2=((1-theta)*(1-fs1)+theta*(1-fh1))*H2;

rho=exp(gam)/disc-(1-delta);

w1=[alpha^alpha*(1-alpha)^(1-alpha)/(1+lambdap)/(rho^alpha)]^(1/(1-alpha));
w2=A2*w1;

klr1=w1/rho*alpha/(1-alpha);
klr2=A2*klr1;
k=klr1*L1;

F1=lambdap/(1+lambdap)*L1*klr1^alpha;
y1=L1*klr1^alpha-F1;
y=y1/v;
y2=(1-v)*y;
F2=lambdap*y2/A2;
x=y;

is=[1-(1-delta)*exp(-gam)]*exp(gam)*k/(1-theta);
g=g_x*x;
c=y-g-(1-theta)*is;

tau=tau_x*x;
tauh1=tau*psih1/theta/fh1;
tauh2=tau*psih2/theta/(1-fh1);
taus=tau-theta*(fh1*tauh1+(1-fh1)*tauh2);

t=t_x*x;
th1=t*psih1/theta/fh1;
th2=t*psih2/theta/(1-fh1);
ts=(t-theta*(fh1*th1+(1-fh1)*th2))/(1-theta);

global w1 w2 H1 H2 th1 th2 tauh1 tauh2 c g t tau
global s theta fh1 fs1 gam piss h disc
ch10=c/theta/fh1;
ch20=c/theta/(1-fh1);
cs0=c/(1-theta);
lambdah10=exp(gam)/(exp(gam)*ch10-h*c);
lambdah20=exp(gam)/(exp(gam)*ch20-h*c);
lambdas0=exp(gam)/(exp(gam)*cs0-h*c);
bR0=(g-t+tau)/(1-(1.025/4)/exp(gam)/piss);
R0=(1.025/4);

% ,'Display','iter'
%options = optimoptions('fsolve','MaxFunctionEvaluations',8000,'MaxIterations',4000,'Display','final-detailed');

ss=fsolve(@distSS,[ch10 ch20 cs0 lambdah10 lambdah20 lambdas0 bR0 R0])
ch1=ss(1);
ch2=ss(2);
cs=ss(3);
lambdah1=ss(4);
lambdah2=ss(5);
lambdas=ss(6);
bR=ss(7);
R=ss(8);


function r=distSS(x)
global w1 w2 H1 H2 th1 th2 tauh1 tauh2 c g t tau
global s theta fh1 fs1 gam piss h disc

ch1=x(1);
ch2=x(2);
cs=x(3);
lambdah1=x(4);
lambdah2=x(5);
lambdas=x(6);
bR=x(7);
R=x(8);

r(1)=w1*H1-th1+tauh1+R*bR*fs1*(1-s)/theta/fh1/exp(gam)/piss-ch1;
r(2)=w2*H2-th2+tauh2+R*bR*(1-fs1)*(1-s)/theta/(1-fh1)/exp(gam)/piss-ch2;
r(3)=c/(1-theta)-theta*(fh1*ch1+(1-fh1)*ch2)/(1-theta)-cs;

r(4)=exp(gam)/(exp(gam)-h)/ch1-lambdah1;
r(5)=exp(gam)/(exp(gam)-h)/ch2-lambdah2;
r(6)=exp(gam)/(exp(gam)-h)/cs-lambdas;
% r(4)=exp(gam)/(exp(gam)*ch1-h*ch1)-lambdah1;  % habit depends on c
% r(5)=exp(gam)/(exp(gam)*ch2-h*ch2)-lambdah2;
% r(6)=exp(gam)/(exp(gam)*cs-h*cs)-lambdas;

r(7)=(g-t+tau)/(1-R/exp(gam)/piss)-bR;

r(8)=exp(gam)*piss/disc*lambdas/[s*lambdas+(1-s)*(fs1*lambdah1+(1-fs1)*lambdah2)]-R;

end
