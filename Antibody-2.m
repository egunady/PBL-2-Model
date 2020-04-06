% Antibody.m

% Calculate plasma blood cell levels
t = linspace(1,28);
Bp = t;
for i = 1:length(t)
    if(t(i) > 6)
        Bp(i) = 20000*(1 - exp(-0.693*(t(i)-6))); %cell/mL

    else 
        Bp = 0;

    end
end

figure(1); clf;
plot(t, Bp);
title('Increase in plasma B cells in the Blood');
xlabel('Time (days)');
ylabel('Concentration (cells/mL)');

%Initial concentrations and constants
times = [0 24];
A0 = 0;
k1 = 200e-12*20000;
k2 = 0.105;


%Calculation of Antibody concentration
[t,A] = ode45(@(t,A) Antibodies(t,A,k1,k2),times, A0);

t= [linspace(0, 6)'; (t+6)]; 
A = [zeros(100,1); A];
threshold = ones(1,length(t))*31.94e-6;

%Plot concentration with average level
figure(2);clf;
plot(t, A, t, threshold, '--')
title('Increase in HIV-1 Antibodies over Acute Phase')
xlabel('Time (days)');
ylabel('Concentration (g/mL)');
legend('[Ig]','Average Concentration');

function dA = Antibodies(t,A,k1,k2)
    %A is X(1), B is X(2), C is X(3), D is X(4)
    dA = 0;
    dA = k1*(1-exp(-0.693*t)) - k2*A; 
end