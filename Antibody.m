% Antibody.m
% Written by: Cooper Lair
% BME 260
% April 6, 2020

% Calculate plasma blood cell levels
t = linspace(1,28);
Bp = t;
for i = 1:length(t)
    if(t(i) > 6)
        Bp(i) = 20*(1 - exp(-0.693*(t(i)-6))); %cell/uL

    else 
        Bp = 0;

    end
end

%Plot B cell concentrations
figure(1); clf;
plot(t, Bp,'color', [0.4940, 0.1840, 0.5560]);
title('Increase in plasma B cells in the Blood');
xlabel('Time (days)');
ylabel('Concentration (molecules/uL)');

%Initial concentrations and constants
times = [0 24];
A0 = 0;
k1 = 200e-12*20; %g/uL
k1 = k1/150000; %Divide my molar mass for mol/mL
k1= 6.022e23*k1; %molecules/uL
k2 = 0.125;


%Calculation of Antibody concentration
[t,A] = ode45(@(t,A) Antibodies(t,A,k1,k2),times, A0);

%Start at Day 6 to account for initial B cell activation
t= [linspace(0, 6)'; (t+6)]; 
A = [zeros(100,1); A];

%Plot line representing average antibody concentration
threshold = ones(1,length(t))*31.94e-9*6.022e23/(150000);

%Plot concentration with average level
figure(2);clf;
plot(t, A, 'm', t, threshold, '--')
title('Increase in HIV-1 Antibodies over Acute Phase')
xlabel('Time (days)');
ylabel('Concentration (molecules/uL)');
legend('[Ig]','Average HIV Positive Level');

function dA = Antibodies(t,A,k1,k2)
    dA = 0;
    dA = k1*(1-exp(-0.693*t)) - k2*A; 
end