function [generated_children,f_best] = dz_2()


k_niz = [2,10,50];
b = 8;    
N = 100; 
max_iter = 20;
f_niz = zeros(3,N+1,max_iter);
ind_best = zeros(1,N);
f_best = zeros(1,N);

for cnt = 1:3 
    k = k_niz(cnt);
    for i = 1:N
        % pocetna random generisanja resenja
        x = -1 + 2*rand(k,1);
        y = -1 + 2*rand(k,1);
        % i odgovarajuce ciljne funkcije
        f = 2-0.1*(cos(5*pi.*x)+cos(5*pi.*y))-x.^2-y.^2;
        % inicijalizacija potrebnih nizova i promenljivih 
        f_new = zeros(1, k*b);
        x_new = zeros(1, k*b);
        y_new = zeros(1, k*b);
        f_global = -100;
        % brojac populacija
        t = 0;

        while t<=max_iter
            t = t+1;
            for j = 1:k
                x_new(1,((j-1)*b+1):b*j) = x(j) + (-0.5 + rand(1,b));
                y_new(1,((j-1)*b+1):b*j) = y(j) + (-0.5 + rand(1,b));
            end
            x_new(x_new>1) = 1; 
            x_new(x_new<-1) = -1; 
            y_new(y_new>1) = 1; 
            y_new(y_new<-1) = -1; 

            f_new = 2-0.1*(cos(5*pi.*x_new)+cos(5*pi.*y_new))-x_new.^2-y_new.^2;
            [f,index] = maxk(f_new,k);
            % najbolja deca --> roditelji u sled krugu
            x = x_new(index);
            y = y_new(index);
            f_niz(cnt,i,t) = mean(f);
            
            % pamcenje najboljeg resenja
            if (f(1) > f_global && cnt == 2)
                f_global = f(1);
                ind_best(i) = t; 
                f_best(i) = f_global;
            end
        end

    end
end

for j = 1:3
    for i = 1:max_iter
        f_niz(j,N+1,i) = mean(f_niz(j,1:N,i));
    end
end
x_max = zeros(3,1);
for j = 1:3
    x_max(j) = nnz(f_niz(j,N+1,:));
end
x_max_val = max(x_max);
%%
figure
for j = 3:-1:1
    fv = f_niz(j,N+1,1:x_max_val);
    fv = reshape(fv,1,x_max_val);
    plot(linspace(1,x_max_val,x_max_val),fv,'LineWidth',1); hold on; 
end
hold off; 
legend('previse','adekvatno','premalo','Location','SouthEast')
xlim([1, x_max_val]); %ylim([1 2.15])
title("ciljna funkcija ~ iteracije")
xlabel("iter"); 

%% racunanje povratnih parametara za histogram

% po generaciji 8 roditelja generise 10 dece
generated_children = k_niz(2)*b*ind_best;


[max(f_niz(1,N+1,:)),max(f_niz(2,N+1,:)),max(f_niz(3,N+1,:))]


end