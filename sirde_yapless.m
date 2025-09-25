clear
close all
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

p_names=['$\beta$','$w$','$\phi$','$\sigma$','$\psi$','$\epsilon$','$\delta$','$S_0$','$I_0$','$R_0$','$D_0$','$E_0$'];

params=[ beta-null, w-null, phi-null, sigma-null, psi-null, epsilon-null, delta-null, S0, I0, R0, D0, E0];
tspan = linspace(0, t-final, intervals);
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