function [f_niz,N_iter,ind_best,f_best] = simulirano_kaljenje(factor)

% ocekivani maks oko 2.13
N = 100; 

if (factor >= 1.1 && factor < 2)
    N_iter = 1000;
elseif (factor < 1.1)
    N_iter = 10000;
else 
    N_iter = 100;
end
f_niz = zeros(N+1,N_iter);
ind_best = zeros(1,N);
f_best = zeros(1,N);

for i = 1:N
    
    k = 0; 
    Tk = 80; 
    f_global_max = -100;   
    
    % generisanje pocetne tacke
    x = -1 + 2*rand();
    y = -1 + 2*rand();
    f = 2-0.1*(cos(5*pi*x)+cos(5*pi*y))-x^2-y^2;
    
    for j = 1:N_iter
        % raspored hladjenja:
        Tk = Tk/factor;
        
        if k<10
            Mk = 1;
        else
            Mk = 5;
        end
        f_curr_max = -100;

        
        for m = 0:Mk
            % generisi x_new,y_new na osn x, u opsegu -1,1
            x_new = x + (-0.5 + rand());
            if (x_new>1)
                x_new = 1;
            elseif (x_new<-1)
                x_new=-1;
            end

            y_new = y + (-0.5 + rand());
            if (y_new>1)
                y_new = 1;
            elseif (y_new<-1)
                y_new=-1;
            end

            f_curr = 2-0.1*(cos(5*pi*x_new)+cos(5*pi*y_new))-x_new^2-y_new^2;
            delta = f_curr - f;
            change = 0; 
            
            if delta>0 
                x = x_new;
                change = 1;
            else
                % za negativno delta
                if exp(delta/Tk) > rand()
                    x = x_new;
                    change = 1;
                end
            end
            if change == 1
                f = f_curr;
                f_curr_max = max(f_curr_max,f);
            end
            
            f_niz(i,j) = f;
        end
        
        if f_curr_max>f_global_max 
            f_global_max = max(f_global_max,f_curr_max);
            ind_best(i) = j;
            f_best(i) = f_global_max;
        end

        k = k + 1;
        
    end

end

% PROMENLJIVE MAX GLOBAL I OVA DRUGA

% title("Vremenska zavisnost uprosecene ciljne funkcije resenja kroz iteracije")
% xlabel("iter"); 