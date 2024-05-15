figure(3)

Vnoise_comparator_matrix_8bits(1)=0;
Vnoise_comparator_matrix_10bits(1)=0;
Vnoise_comparator_matrix_12bits(1)=0;

errorbar(Vnoise_comparator_matrix_8bits,INL_Vnoise_Comparator_8bits,INL_Vnoise_Comparator_tolerance_8bits,'o','LineWidth',1.01)
hold on
errorbar(Vnoise_comparator_matrix_10bits,INL_Vnoise_Comparator_10bits,INL_Vnoise_Comparator_tolerance_10bits,'o','LineWidth',1.01)
errorbar(Vnoise_comparator_matrix_12bits,INL_Vnoise_Comparator_12bits,INL_Vnoise_Comparator_tolerance_12bits,'o','LineWidth',1.01)

plot(Vnoise_comparator_matrix_8bits,INL_Vnoise_Comparator_8bits,'-k','LineWidth',1.01)
plot(Vnoise_comparator_matrix_10bits,INL_Vnoise_Comparator_10bits,'--k','LineWidth',1.01)
plot(Vnoise_comparator_matrix_12bits,INL_Vnoise_Comparator_12bits,'-.k','LineWidth',1.01)
set(gca, 'XScale','log')
% for k = [1 2 3] 
%     errorbar(Vnoise_comparator_matrix,INL_Vnoise_Comparator_together(k,:),INL_Vnoise_Comparator_tolerance_together(k,:),allmarks{k},'LineWidth',1.01);
%     hold on
% end

title('INL_{máx} vs Vnoise\_comparator')
xlabel('Vnoise\_comparator [V]')
ylabel('INL_{máx} [LSB]')
legend('8 bits - Variable T-Step','10 bits - Variable T-Step','12bits - Variable T-Step','8 bits - Fixed T-Step','10 bits - Fixed T-Step','12bits - Fixed T-Step')
grid on
hold off

%%

for k = [1 2 3]
    errorbar(e2_matrix*100,INL_e2_together(k,:),INL_e2_tolerance_together(k,:),'o','LineWidth',1.01);
    hold on
end

plot(e2_matrix*100,INL_e2_together(1,:),'-k','LineWidth',1.01);
plot(e2_matrix*100,INL_e2_together(2,:),'--k','LineWidth',1.01);
plot(e2_matrix*100,INL_e2_together(3,:),'-.k','LineWidth',1.01);
grid on
title('INL_{máx} vs e2')
xlabel('e2 [%/V^2]')
ylabel('INL_{máx} [LSB]')
legend('8 bits - Variable T-Step','10 bits - Variable T-Step','12bits - Variable T-Step','8 bits - Fixed T-Step','10 bits - Fixed T-Step','12bits - Fixed T-Step')
hold off

%%

figure(14)
ENOB_tau_together(3,1:1:3)=11.84;
for k = [1 2 3]
    plot(tau_matrix,ENOB_tau_together(k,:),'o','LineWidth',1.01);    
    hold on
end

plot(tau_matrix,ENOB_tau_together(1,:),'-k','LineWidth',1.01);  
plot(tau_matrix,ENOB_tau_together(2,:),'--k','LineWidth',1.01);  
plot(tau_matrix,ENOB_tau_together(3,:),'-.k','LineWidth',1.01);  


grid on
title('ENOB vs \tau_{comp}')
xlabel('\tau_{comp} [s]')
ylabel('ENOB [bits]')
legend('8 bits - Variable T-Step','10 bits - Variable T-Step','12bits - Variable T-Step','8 bits - Fixed T-Step','10 bits - Fixed T-Step','12bits - Fixed T-Step')

hold off

%%

%% JITTER CON RELOJ EN TIEMPO CONTINUO - V4 (+PUNTOS, + DIVERSIÓN)- LIMITED EDITION

allmarks = {'o-','*-','x-','+-','.-','s-','d-','^-','v-','>-','<-','p-','h-'};
% GHZ= 0.98* (1/10)*([0.1,1,10,100,1000,10000]*10^3)/(2^B); % represento fin
% 20 ns
% ps2_ENOB = [10,10,10,9.94, 9.12,4.50];
% ps200_ENOB = [10, 10, 9.12, 4.50];
% ps5_ENOB = [10,10,10,9.88, 8.05, 3.58];
% ps10_ENOB = [];
% ps20_ENOB = [10,10, 10, 9.12, 4.5];
% ps50_ENOB = [10,10, 9.88, 8.05, 3.58];
% ps100_ENOB = [10,10, 9.66, 5.99, 2.91];

int1 = [0.1,0.2,0.4,0.6,0.8];
int2 = [1,2,4,6,8];
int3 = [10,20,40,60,80];
int4 = [100,200,400,600,800];
int5 = [1000,2000,4000,6000,8000];
int6 = [10000,20000];

GHZ = 0.98 * (1/10)*([int1,int2,int3,int4,int5,int6]*10^3)/(2^10);
%
ps5_int1 = [9.97,9.97,9.97,9.97,9.97];
ps5_int2 = [9.97,9.97,9.97,9.97,9.97];
ps5_int3 = [9.97,9.97,9.96,9.94,9.92];
ps5_int4 = [9.88,9.69,9.17,8.73,8.36];
ps5_int5 = [8.09,6.23,4.72,4.27,3.85];
ps5_int6 = [3.56,2.79];
%
ps10_int1 = [9.97,9.97,9.97,9.97,9.97];
ps10_int2 = [9.97,9.97,9.97,9.97,9.97];
ps10_int3 = [9.96,9.96,9.92,9.85,9.79];
ps10_int4 = [9.69,9.17,8.36,6.94,6.67];
ps10_int5 = [6.23,4.72,3.85,3.37,3.09];
%
ps20_int1 = [9.97,9.97,9.97,9.97,9.97];
ps20_int2 = [9.97,9.97,9.97,9.96,9.96]; 
ps20_int3 = [9.96,9.92,9.79,9.59,9.39];
ps20_int4 = [9.17,8.36,6.67,5.39,4.95];
ps20_int5 = [4.72,3.85,3.09];
%
ps50_int1 = [9.97,9.97,9.97,9.97,9.97];
ps50_int2 = [9.97,9.97,9.96,9.94,9.92];
ps50_int3 = [9.88,9.69,9.17,8.73,8.36];
ps50_int4 = [8.09,6.23,4.72,4.27,3.85];
ps50_int5 = [3.56,2.79];
%
ps100_int1 = [9.97,9.97,9.97,9.97,9.97];
ps100_int2 = [9.96,9.96,9.92,9.85,9.79];
ps100_int3 = [9.69,9.17,8.36,6.94,6.67];
ps100_int4 = [6.23,4.72,3.85,3.37,3.09];

ps5_ENOB = [ps5_int1,ps5_int2,ps5_int3,ps5_int4,ps5_int5,ps5_int6];
ps10_ENOB = [ps10_int1,ps10_int2,ps10_int3,ps10_int4,ps10_int5];
ps20_ENOB = [ps20_int1,ps20_int2,ps20_int3,ps20_int4,ps20_int5];
ps50_ENOB = [ps50_int1,ps50_int2,ps50_int3,ps50_int4,ps50_int5];
ps100_ENOB = [ps100_int1,ps100_int2,ps100_int3,ps100_int4];




plot(GHZ(1:23),ps5_ENOB(1:23),'o','LineWidth',1.01)
hold on
plot(GHZ(1:22),ps10_ENOB(1:22),'o','LineWidth',1.01)
plot(GHZ(1:18),ps50_ENOB(1:18),'o','LineWidth',1.01)
plot(GHZ(1:17),ps100_ENOB(1:17),'o','LineWidth',1.01)

plot(GHZ(1:23),ps5_ENOB(1:23),'-k','LineWidth',1.01)
hold on
plot(GHZ(1:22),ps10_ENOB(1:22),'--k','LineWidth',1.01)
plot(GHZ(1:18),ps50_ENOB(1:18),'.-k','LineWidth',1.01)
plot(GHZ(1:17),ps100_ENOB(1:17),':k','LineWidth',1.01)



% plot(GHZ(1:4),ps200_ENOB,allmarks{5},'LineWidth',1.01)

set(gca, 'XScale','log')

title('ENOB for differente RMS Aperture Jitter Values - ADC 10 bits')
xlabel('Analog input frequency [MHz]')
ylabel('ENOB [bits]')
legend('5ps - Variable T-Step','10 ps - Variable T-Step','50 ps - Variable T-Step','100 ps - Variable T-Step','5ps - Fixed T-Step','10 ps - Fixed T-Step','50 ps - Fixed T-Step','100 ps - Fixed T-Step')
grid on
hold off



