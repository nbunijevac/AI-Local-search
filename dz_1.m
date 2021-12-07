function [generated_children,f_best] = dz_1()
 
N = 100;

%% dobro hladjenje 

[f_niz1,N_iter,ind_best,f_best] = simulirano_kaljenje(1.1);

for i = 1:N_iter
    f_niz1(N+1,i) = mean(f_niz1(1:N,i));
end

x_max = nnz(f_niz1(N+1,:));
figure(1)
plot(linspace(1,N_iter,N_iter),f_niz1(N+1,:),'LineWidth',1)
xlim([1, x_max]); ylim([1 2.2])
title("Vremenska zavisnost uprosecene ciljne funkcije resenja kroz iteracije")
legend('adekvatno hladjenje', 'Location', 'East')
xlabel("iter"); 

%% brze hladnjenje

[f_niz2,N_iter,~,~] = simulirano_kaljenje(200);

for i = 1:N_iter
    f_niz2(N+1,i) = mean(f_niz2(1:N,i));
end

x_max = nnz(f_niz2(N+1,:));
figure(2)
plot(linspace(1,N_iter,N_iter),f_niz2(N+1,:),'LineWidth',1)
xlim([1, x_max]); ylim([1 2.15])
title("Brze hladjenje")
xlabel("iter"); 

%% sporije hladjenje

[f_niz3,N_iter,~,~] = simulirano_kaljenje(1.001);

for i = 1:N_iter
    f_niz3(N+1,i) = mean(f_niz3(1:N,i));
end

x_max = nnz(f_niz3(N+1,:));
figure(3)
plot(linspace(1,N_iter,N_iter),f_niz3(N+1,:),'LineWidth',1)
xlim([1, x_max]); ylim([1 2.15])
title("Sporije hladjenje")
xlabel("iter"); 




%% modifikacija

[f_niz,N_iter,~] = simulirano_kaljenje_modifikacija(1.1);

for i = 1:N_iter
    f_niz(N+1,i) = mean(f_niz(1:N,i));
end

x_max = nnz(f_niz(N+1,:));

figure(4)
plot(linspace(1,N_iter,N_iter),f_niz1(N+1,:),'LineWidth',1)
hold on
plot(linspace(1,N_iter,N_iter),f_niz(N+1,:),'LineWidth',1)
xlim([1, x_max]); ylim([1 2.15])
legend('original', 'modifikacija','Location','SouthEast')
title("ciljna funkcija ~ iteracije")
xlabel("iter"); 

[max(f_niz1(N+1,:)),max(f_niz2(N+1,:)),max(f_niz3(N+1,:)),max(f_niz(N+1,:))]


%% izracunavanje povratnih parametara za histogram

% po jedno dete po generaciji
generated_children = ind_best;



end
