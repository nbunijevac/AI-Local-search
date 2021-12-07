clc
clear all
close all

[generated_children_sk,f_sk] = dz_1();
[generated_children_pps,f_pps] = dz_2();
[generated_children_ga,f_ga] = dz_3();
N_stub = 20;

figure
subplot(3,1,1), histogram(generated_children_sk,N_stub)
title('Broj generisanih cvorova do najboljeg resenja za:')
ylabel('SK')
subplot(3,1,2), histogram(generated_children_pps,N_stub)
ylabel('PPS')
subplot(3,1,3), histogram(generated_children_ga,N_stub)
ylabel('GA')


figure
subplot(3,1,1), histogram(f_sk,N_stub)
title('Vrednosti ciljne f-je u trenutku najboljeg resenja za:')
ylabel('SK')
subplot(3,1,2), histogram(f_pps,N_stub)
ylabel('PPS')
subplot(3,1,3), histogram(f_ga,N_stub)
ylabel('GA')
