function [generated_children,f_best] = dz_3()

k_niz = [2,50,200];
    
N = 100; 
max_iter = 200;
f_niz = zeros(3,N+1,max_iter);
p_mutation = 0.001;
ind_best = zeros(1,N);
f_best = zeros(1,N);


for cnt = 1:3 
    k = k_niz(cnt);
    for i = 1:N
        f_global = -100;
        
        % pocetna random generisanja resenja
        x = -1 + 2*rand(1,k);
        y = -1 + 2*rand(1,k);
        % i odgovarajuce ciljne funkcije
        f = 2-0.1*(cos(5*pi.*x)+cos(5*pi.*y))-x.^2-y.^2;
        [f,index] = sort(f,'descend');
        % kodiracemo koordinate kao pozitivne 8-bitne int vrednosti
        % (-1,0) -> (0,5e7)
        % (0,1)  -> (5e7+1,1e8-1) 
        x = floor(x(index)*5e7+5e7);
        y = floor(y(index)*5e7+5e7);
        x(x==1e8) = 1e8-1;
        y(y==1e8) = 1e8-1;
        % niz procenata
        f_proc = f/sum(f);
        f_proc = cumsum(f_proc);
        % inicijalizacija potrebnih nizova i promenljivih 
        x_new = zeros(1,k);
        y_new = zeros(1,k);
        % brojac populacija
        t = 0;
        
        while t <= max_iter
            t = t+1;
            % kompletno generisanje potomaka 
            for j = 1:k
                % prvi roditelj
                num = rand();
                index_f_curr = k - nnz((f_proc-num)>0) + 1;
                first_parent_x = x(index_f_curr);
                first_parent_y = y(index_f_curr);
                % drugi roditelj
                num = rand();
                index_f_curr = k - nnz((f_proc-num)>0) + 1;
                second_parent_x = x(index_f_curr);
                second_parent_y = y(index_f_curr);
                % UKRSTANJE
                % cutoff sa desne strane 
                cutoff = randi(8,1);
                x_new(j) = floor(floor(first_parent_x/10^cutoff)*10^cutoff + abs(rem(second_parent_x,10^cutoff)));
                y_new(j) = floor(floor(first_parent_y/10^cutoff)*10^cutoff + abs(rem(second_parent_y,10^cutoff)));
                % MUTACIJA
                for ind = 0:7
                    if (rand()<p_mutation)
                        x_new(j) =  (floor(x_new(j)/10^(ind+1))*10 + randi(10)-1)*10^(ind) + mod(x_new(j),10^(ind));
                    end
                    if (rand()<p_mutation)
                        y_new(j) =  (floor(y_new(j)/10^(ind+1))*10 + randi(10)-1)*10^(ind) + mod(y_new(j),10^(ind));
                    end
                end
            end
            
            % deca -> roditelji u sled krugu
            x = x_new;
            y = y_new;

            x(x>=1e8) = 1e8-1;
            y(y>=1e8) = 1e8-1;

            % racunanje ciljne f-je
            f = 2-0.1*(cos(5*pi.*((x-5e7)/5e7))+cos(5*pi.*((y-5e7)/5e7)))-((x-5e7)/5e7).^2-((y-5e7)/5e7).^2;
            % sortiranje
            [f,index] = sort(f,'descend');
            x = x(index);
            y = y(index);
            % racunanje procenata
            f_proc = f/sum(f);
            f_proc = cumsum(f_proc);

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
legend('previse','adekvatno','premalo','Location','East')
xlim([0, x_max_val]); %ylim([1 2.15])
title("ciljna funkcija ~ iteracije")
xlabel("iter"); 

%% racunanje povratnih parametara za histogram

% po generaciji 8 roditelja generise 10 dece
generated_children = k*ind_best;

[max(f_niz(1,N+1,:)),max(f_niz(2,N+1,:)),max(f_niz(3,N+1,:))]


end