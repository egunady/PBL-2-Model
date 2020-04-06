%updated code with lymph and blood transport rates + cd8, cytokines
tspan = [0 35];
%HIV Infected Initial Conditions
%IC = [1.05E+10 0 5.25 5.0E+9 0 0 4000 0];
%Healthy Initial Conditions
IC = [1.05E+10 0 0 5.0E+9 0 0 0 0];

%Colors
orange = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];
set(0,'defaultAxesFontSize',14)

[t, c] = ode23s(@conc, tspan, IC);
figure(1);
clf();
yyaxis left;
plot(t, (1/1.05E+07)*c(:,1), t, (1/1.05E+07)*c(:,2), 'LineWidth', 1);
ylim([0 1500]);
ylabel('Concentration of T cells (cells/uL)');
yyaxis right;
plot(t, (1/1.05E+07)*c(:,3),'LineWidth', 1);
ylabel(' Concentration of HIV (virions/uL)');
xlabel('Time, days');
xlim([0 30]);
ylim([0 1]);
legend('Healthy CD4+ T cells', 'Infected CD4+ T cells', 'HIV');
title('Concentration of CD4+ T cells and HIV in lymphatic tissue');

figure(2);
clf();
yyaxis left;
plot(t, (1/5E+06)*c(:,4), t, (1/5E+06)*c(:,5), 'LineWidth', 1);
ylim([0 2000]);
ylabel('Concentration of T cells (cells/uL)');
yyaxis right;
plot(t, (1/5E+06)*c(:,6),'LineWidth', 1);
ylabel('Concentration of HIV (virions/uL)');
xlabel('Time, days');
xlim([0 30]);
ylim([0 1]);
legend('Healthy CD4+ T cells', 'Infected CD4+ T cells', 'HIV');
title('Concentration of CD4+ T cells and HIV in blood');

fig3 = figure(3);
clf();
set(fig3, 'defaultAxesColorOrder', [orange; green]);
yyaxis left;
plot(t, (1/1.05E+07)*c(:,3), 'Color', orange, 'LineWidth', 1);
ylabel('Concentraiton of HIV, (virions/uL)', 'Color', orange);
yyaxis right;
plot(t, (1/1.05E+07)*c(:,7), 'Color', green, 'LineStyle', '--', 'LineWidth', 1);
ylabel('Concentration of Il-1beta (proteins/uL)');
xlim([5 15]);
legend('HIV', 'Il-1\beta');
xlabel('Time, days');
title('Concentration of HIV and Il-1\beta in lymphatic tissue');

fig4 = figure(4);
clf();
set(fig4, 'defaultAxesColorOrder', [0 0 0; orange]);
yyaxis left;
plot(t, (1/5E+06)*c(:,6), 'k-', t, (1/1.05E+07)*c(:,3), 'k--','LineWidth', 1)
ylim([0 8E+05]);
ylabel('Concentration of HIV (virions/uL)');
yyaxis right;
plot(t, (1/10500)*c(:,8), '-.', 'Color', orange, 'LineWidth', 1);
xlim([5 15]);
ylabel('Concentration of HIV (virions/uL)');
xlabel('Time, days');
title('Concentration of HIV');
legend('Blood', 'Paracortex', 'Follicular Dendritic Cells');


function [dcdt] = conc(t, c)
dcdt = zeros(7, 1);
dT = 0.01;
r1 = 0.0002;
r2 = 0.01;
lambda1 = dT*5E+09 - r1*5E+09 + r2*1E+08;
%lambda2 = 0; %HIV Infected
lambda2 = 1E+08; %Healthy
m1 = 0.1;
m2 = 0.2;
beta = 2.29E-12;
delta = 1;
p = 2.5E+4;
clear = 100;
k1 = 200e-12*20000;
k2 = 0.105; 
Nc = 15;
d3 =0.001;
d5 = 6.6;
%Paracortex: c(1) = healthy T cells, c(2) = infected T cells, c(3) = viral load
%Blood: c(4) = healthy T cells, c(5) = infected T cells, c(6) = viral load
%c(7) = cytokines in lymph (doesn't accumulate in blood)
%c(8) = HIV in follicular dendritic cells
dcdt(1) = lambda1 - dT*c(1) - beta*c(3)*c(1) - r1*c(1) + r2*c(4); 
dcdt(2) = beta*c(3)*c(1) - delta*c(2) - r1*c(2) + r2*c(5);
dcdt(3) = p*c(2) - clear*c(3) - m1*c(3) + m2*c(6);
dcdt(4) = lambda2 - dT*c(4) - beta*c(6)*c(4) - r2*c(4) + r1*c(1); 
dcdt(5) = beta*c(6)*c(4) - delta*c(5) - r2*c(5) + r1*c(2);
dcdt(6) = p*c(5) - clear*c(6) - m2*c(6) + m1*c(3);
dcdt(7) = Nc*d3*(c(2)) - d5*c(7);
dcdt(8) = c(3) - m1*c(3);

end