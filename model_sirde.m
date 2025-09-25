%Ebola Model
%untested & dataless </3

clear
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

p_names=['$\beta$','$w$','$\phi$','$\sigma$','$\psi$','$\epsilon$','$\delta$','$S_0$','$I_0$','$R_0$','$D_0$','$E_0$'];

%We may change the parameter set to observe behavioral changes of Ebola in our population
params=[ beta-null, w-null, phi-null, sigma-null, psi-null, epsilon-null, delta-null, S0, I0, R0, D0, E0];
%Define right hand ODE solution set on some time interval with the 5 initial conditions
tspan = linspace(0, t-final, intervals);
%t-final will be an estimation of when our endemic "ends" which we're able to figure out by increasing t arbitrarily to test steady states (limits) or grabbing average endemic lengths.
[t,y_solution]=ode45(@(t,Y) rhs(t,Y,params), [tspan(1),tspan(end)], [params(end-4),params(end-3),params(end-2),params(end-1),params(end)]); 
figure(1)

plot(t,y_solution(:,1), 'LineWidth',2)
hold on
plot(t,y_solution(:,2), 'LineWidth',2)
plot(t,y_solution(:,3), 'LineWidth',2)
plot(t,y_solution(:,4), 'LineWidth',2)
plot(t,y_solution(:,5), 'LineWidth',2)

legend(['Susceptible','Infected','Recovered','Deceased','Exposed'], fontsize=16)
xlabel('Time', fontsize=16)
ylabel('Population', fontsize=16)

function f1 = rhs_S(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    D=Y(4);
    E=Y(5);
    beta=params(1);
    w=params(2);
    phi=params(3);
    sigma=params(4);
    psi=params(5);
    epsilon=params(6);
    delta=params(7);
    S0=params(8);
    I0=params(9);
    R0=params(10);
    D0=params(11);
    E0=params(12);
    f1 = -sigma*D*S-beta*S*E + phi*R;
end

function f2 = rhs_I(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    D=Y(4);
    E=Y(5);
    beta=params(1);
    w=params(2);
    phi=params(3);
    sigma=params(4);
    psi=params(5);
    epsilon=params(6);
    delta=params(7);
    S0=params(8);
    I0=params(9);
    R0=params(10);
    D0=params(11);
    E0=params(12);
    f2 = epsilon*E-delta*I-w*I;
end

%Note, I added a R->S parameter as phi, this appears symmetrically through f1
%phi=0 in our parameter space defined above removes this variable
function f3 = rhs_R(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    D=Y(4);
    E=Y(5);
    beta=params(1);
    w=params(2);
    phi=params(3);
    sigma=params(4);
    psi=params(5);
    epsilon=params(6);
    delta=params(7);
    S0=params(8);
    I0=params(9);
    R0=params(10);
    D0=params(11);
    E0=params(12);
    f3 = w*I-phi*R;
end

function f4 = rhs_D(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    D=Y(4);
    E=Y(5);
    beta=params(1);
    w=params(2);
    phi=params(3);
    sigma=params(4);
    psi=params(5);
    epsilon=params(6);
    delta=params(7);
    S0=params(8);
    I0=params(9);
    R0=params(10);
    D0=params(11);
    E0=params(12);
    f4 = psi*E+delta*I;
end

%so the rhsY(n) relates to some categorical variable (Susceptible,Dead,etc)
function f5 = rhs_E(t,Y,params)
    S=Y(1);
    I=Y(2);
    R=Y(3);
    D=Y(4);
    E=Y(5);
    beta=params(1);
    w=params(2);
    phi=params(3);
    sigma=params(4);
    psi=params(5);
    epsilon=params(6);
    delta=params(7);
    S0=params(8);
    I0=params(9);
    R0=params(10);
    D0=params(11);
    E0=params(12);
    f5 = sigma*D*S+beta*S*E-psi*E-epsilon*I;
end

%Suppose we want probability infected recover, and suppose infected get
%treated with some medicine. Then, using study 
%pubmed.ncbi.nlm.nih.gov/31774950/ treatment, REGN-EB3 mortality
%rate at 4 weeks is .335, we can estimate probability(infected_treated=dead)=.335
%and probability(infected_treated=recover)=.665. Alright, with 1975 kikwit zaire data,
%about 81% died untreated so probability(infected_untreated=dead)=.81. These
%will be interesting numbers that pushes us to test parameters with these
%in mind. That difference of mortality between treated and untreated tells 
%us the death population is significantly different in either case.

%Defining Exposed to contain those Infected_untreated allows us to
%use the treatment probabilities for the Infected, as all Infected become
%Infected_Treated (change of variables). This should carry throughout the
%model (affects relations between categories) but it'll allow us to use data

%How does ^this^ affect healthcare workers? Well, E and I respectively become
%E = {ill and untreated} and I = {ill and treated}. This also allows us to 
%remove the healthcare*exposed as those untreated won't come in contact by definition
%This allows us to predict if the virus may become endemic with respect to
%medicine

%What'd I add? Recover to Susceptible. What do I suggest? E&I redefinition.