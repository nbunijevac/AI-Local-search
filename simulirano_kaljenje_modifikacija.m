function [f_niz,N_iter,ind_best] = simulirano_kaljenje_modifikacija(factor)

N = 100; 

if (factor == 1.1)
    N_iter = 1000;
elseif (factor < 1.1)
    N_iter = 10000;
else 
    N_iter = 100;
end
f_niz = zeros(N+1,N_iter);
f_global_max = -100;   
k_best = 1;
x_best = -2;
y_best = -2;
ind_best = zeros(1,N);
f_best = zeros(1,N);

for i = 1:N
    
    Tk = 80; 
    
    x = -1 + 2*rand();
    y = -1 + 2*rand();
    f = 2-0.1*(cos(5*pi*x)+cos(5*pi*y))-x^2-y^2;
    
    for k = 1:N_iter
        % raspored hladjenja:
        Tk = Tk/factor;
        
        if k<floor(N_iter/20)
            Mk = 1;
        else
            Mk = 5;
        end

        if (k - k_best > floor(N_iter/50))
            x = x_best;
            y = y_best;
        end
        
        % max na curr temperaturi 
        f_curr_max = -100;
        for m = 0:Mk
            % generisi x_new,y_new na osn x, u opsegu -1,1
            x_new = x + (-0.5 + rand());
            if (x_new>1)
                x_new = 1;
            elseif (x_new<-1)
                x_new = -1;
            end

            y_new = y + (-0.5 + rand());
            if (y_new>1)
                y_new = 1;
            elseif (y_new<-1)
                y_new = -1;
            end

            f_curr = 2-0.1*(cos(5*pi*x_new)+cos(5*pi*y_new))-x_new^2-y_new^2;
            delta = f_curr - f;
            change = 0; 
            
            if delta>0 
                x = x_new;
                change = 1;
            else
                % za negativno delta
                if exp(delta/Tk)> rand()
                    x = x_new;
                    change = 1;
                end
            end
            if change == 1
                f = f_curr;
                f_curr_max = max(f_curr_max,f);
                
            end
            
            f_niz(i,k) = f;
        end
        
        if (f_curr_max > f_global_max)
            k_best = k;
            x_best = x;
            y_best = y;
            f_global_max = f_curr_max;
            ind_best(i) = k;
            f_best(i) = f_global_max;
        end 
        
    end

end