%Replication of Cooper (1999)
%Saeed Cherik
%ECO 6363
clear; clc;


beta= 0.9;
lambda= 0.75;
F= 0.2;
delta = 0.1;

AL= 0.75;
AH= 1.25;
A_grid= [AL, AH];
NA= length(A_grid);

P_A= [0.9 0.1;
       0.1 0.9];

eps_min= 0.4;
eps_max= 1.6;
Neps= 20;
eps_grid= linspace(eps_min,eps_max,Neps);

Nk= 50;
k_grid= (1 - delta).^(0:(Nk-1));
k_grid= k_grid';

revenue= zeros(Nk,NA,Neps);



V= zeros(Nk,NA,Neps);

VR= zeros(Nk,NA,Neps);
VNR= zeros(Nk,NA,Neps);
z= zeros(Nk,NA,Neps);

for ik= 1:Nk
    for ia= 1:NA
        for ie= 1:Neps
            revenue(ik,ia,ie)= A_grid(ia) * eps_grid(ie) * k_grid(ik);
        end
    end
end

mu= 1e-6;
max_iter= 1000;
diff= 1;
iter= 0;

while diff > mu && iter < max_iter
    iter= iter + 1;
    
    EV_replace= zeros(NA,1);
    
    for ia= 1:NA
        temp_sum= 0;
        for ja= 1:NA
           
            avg_over_eps= 0;
            for je= 1:Neps
                avg_over_eps= avg_over_eps + V(1,ja,je);
            end
            avg_over_eps= avg_over_eps / Neps;
            
            
            temp_sum= temp_sum + P_A(ia,ja) * avg_over_eps;
        end
        EV_replace(ia)= temp_sum;
    end
    

    EV_norepl= zeros(Nk,NA);
    
    for ik= 1:Nk

        if ik < Nk
            ik_next= ik + 1;
        else
            ik_next= Nk;
        end
        
        for ia= 1:NA
            temp_sum= 0;
            for ja= 1:NA
                
                avg_over_eps= 0;
                for je= 1:Neps
                    avg_over_eps= avg_over_eps + V(ik_next,ja,je);
                end
                avg_over_eps= avg_over_eps / Neps;
                
                temp_sum= temp_sum + P_A(ia,ja) * avg_over_eps;
            end
            EV_norepl(ik,ia)= temp_sum;
        end
    end

    
    new_V= zeros(Nk,NA,Neps);
    
    for ik= 1:Nk
        for ia= 1:NA
            for ie= 1:Neps
                
                current_rev_NR= A_grid(ia) * eps_grid(ie) * k_grid(ik);
                current_rev_R =A_grid(ia) * eps_grid(ie) * lambda * k_grid(ik) - F;
                
                VR(ik,ia,ie) = current_rev_R  + beta * EV_replace(ia);
                VNR(ik,ia,ie)= current_rev_NR + beta * EV_norepl(ik,ia);
                
                if VR(ik,ia,ie) >= VNR(ik,ia,ie)
                    new_V(ik,ia,ie) = VR(ik,ia,ie);
                    z(ik,ia,ie) = 1;
                else
                    new_V(ik,ia,ie) = VNR(ik,ia,ie);
                    z(ik,ia,ie) = 0;
                end
                
            end
        end
    end
     
    diff= max(abs(new_V(:)-V(:)));
    V= new_V;
    
end

always_replace = zeros(Nk,1);

for ik= 1:Nk
    zz= reshape(z(ik,:,:),NA*Neps,1);
    if all(zz == 1)
        always_replace(ik)= 1;
    end
end

cutoff_index= 0;

for ik= Nk:-1:1
    if always_replace(ik) == 1
        cutoff_index= ik;
    else
        if cutoff_index > 0
            break;
        end
    end
end

if cutoff_index == 0
    disp('No clear cut-off found: there is no k where replacement is always chosen.');
else
    fprintf('Cut-off index for capital (always replace at and below this k): %d\n', cutoff_index);
    fprintf('Cut-off capital level k_cut = %f\n', k_grid(cutoff_index));
end

if cutoff_index > 0
    
    Nk_small= cutoff_index;
    k_grid_small= k_grid(1:Nk_small);

    Nk= Nk_small;
    k_grid= k_grid_small;

    revenue= zeros(Nk,NA,Neps);
    V= zeros(Nk,NA,Neps);
    VR= zeros(Nk,NA,Neps);
    VNR= zeros(Nk,NA,Neps);
    z= zeros(Nk,NA,Neps);

    for ik= 1:Nk
        for ia= 1:NA
            for ie= 1:Neps
                revenue(ik,ia,ie)= A_grid(ia) * eps_grid(ie) * k_grid(ik);
            end
        end
    end

    diff= 1;
    iter= 0;
    max_iter= 1000;

    while diff > mu && iter < max_iter
        iter = iter + 1;

        EV_replace= zeros(NA,1);
        for ia= 1:NA
            temp_sum= 0;
            for ja= 1:NA
                avg_over_eps= 0;
                for je= 1:Neps
                    avg_over_eps = avg_over_eps + V(1,ja,je);
                end
                avg_over_eps= avg_over_eps / Neps;
                temp_sum= temp_sum + P_A(ia,ja) * avg_over_eps;
            end
            EV_replace(ia)= temp_sum;
        end

        EV_norepl= zeros(Nk,NA);
        for ik = 1:Nk
            if ik < Nk
                ik_next= ik + 1;
            else
                ik_next= Nk;
            end
            for ia= 1:NA
                temp_sum= 0;
                for ja= 1:NA
                    avg_over_eps= 0;
                    for je= 1:Neps
                        avg_over_eps = avg_over_eps + V(ik_next,ja,je);
                    end
                    avg_over_eps= avg_over_eps / Neps;
                    temp_sum= temp_sum + P_A(ia,ja) * avg_over_eps;
                end
                EV_norepl(ik,ia)= temp_sum;
            end
        end

        new_V= zeros(Nk,NA,Neps);

        for ik= 1:Nk
            for ia= 1:NA
                for ie= 1:Neps

                    current_rev_NR= A_grid(ia) * eps_grid(ie) * k_grid(ik);
                    current_rev_R= A_grid(ia) * eps_grid(ie) * lambda * k_grid(ik) - F;

                    VR(ik,ia,ie)= current_rev_R  + beta * EV_replace(ia);
                    VNR(ik,ia,ie)= current_rev_NR + beta * EV_norepl(ik,ia);

                    if VR(ik,ia,ie) >= VNR(ik,ia,ie)
                        new_V(ik,ia,ie) = VR(ik,ia,ie);
                        z(ik,ia,ie) = 1;
                    else
                        new_V(ik,ia,ie) = VNR(ik,ia,ie);
                        z(ik,ia,ie) = 0;
                    end

                end
            end
        end

        diff= max(abs(new_V(:) - V(:)));
        V= new_V;
    end

    fprintf('Finished VFI with reduced k-grid: Nk = %d\n', Nk);

end



figure;

subplot(1,2,1);
ia= 1;
z_low= squeeze(z(:,ia, :));

spy(z_low);
title('Policy Function (z=1) for A = A_L');
xlabel('Idiosyncratic productivity index (epsilon)');
ylabel('Capital index (k)');

subplot(1,2,2);
ia= 2;
z_high= squeeze(z(:,ia, :));

spy(z_high);
title('Policy Function (z=1) for A = A_H');
xlabel('Idiosyncratic productivity index (epsilon)');
ylabel('Capital index (k)');


hazard= zeros(Nk,NA);

for ik= 1:Nk
    for ia= 1:NA

        sum_z = 0;
        for ie= 1:Neps
            sum_z = sum_z + z(ik,ia,ie);
        end
        hazard(ik,ia) = sum_z / Neps;
    end
end

age_vec= (0:(Nk-1))';

figure;
plot(age_vec,hazard(:,1), 'o-');
hold on;
plot(age_vec,hazard(:,2), 'x-');
hold off;

xlabel('Capital age (periods since replacement)');
ylabel('Hazard of replacement');
title('Hazard function of capital replacement by aggregate state A');
legend('A = A_L (low)', 'A = A_H (high)', 'Location', 'best');

T= 40;

A_index= zeros(T+1,1);
eps_index= zeros(T+1,1);

k_series= zeros(T+1,1);
A_series= zeros(T+1,1);
 eps_series= zeros(T+1,1);
y_series= zeros(T,1);

k_index(1)= 1;
A_index(1)= 1;
eps_index(1)= ceil(Neps/2);

k_series(1)= k_grid(k_index(1));
A_series(1)= A_grid(A_index(1));
eps_series(1)= eps_grid(eps_index(1));

for t= 1:T
    
    ik= k_index(t);
    ia= A_index(t);
    ie= eps_index(t);
    
    k_t= k_grid(ik);
    A_t= A_grid(ia);
    eps_t= eps_grid(ie);
    

    y_series(t)= A_t * eps_t * k_t;
    
    z_t = z(ik,ia,ie);

    if z_t == 1
        k_index(t+1) = 1;
    else
        
        if ik < Nk
            k_index(t+1) = ik + 1;
        else
            k_index(t+1) = Nk;
        end
    end
    
    k_series(t+1)= k_grid(k_index(t+1));
    
    if t < T+1
        u= rand;
        curA= A_index(t);
        
        if u <= P_A(curA,1)
            A_index(t+1) = 1;
        else
            A_index(t+1) = 2;
        end
        
        A_series(t+1) = A_grid(A_index(t+1));
    end
    
    if t < T+1
        eps_index(t+1) = randi(Neps);
        eps_series(t+1) = eps_grid(eps_index(t+1));
    end
    
end

time_vec = 1:T;

figure;
subplot(2,1,1);
plot(time_vec, k_series(1:T), '-o');
xlabel('Time');
ylabel('Capital k_t');
title('Simulated capital stock for one firm');

subplot(2,1,2);
plot(time_vec, y_series, '-o');
xlabel('Time');
ylabel('Output y_t');
title('Simulated output for one firm');


ia_fixed= 2; 

h= hazard(:, ia_fixed);

dist_k= ones(Nk, 1) / Nk;
T_agg= 40;

I_path= zeros(T_agg, 1);

for t= 1:T_agg
    
    I_path(t) = sum(dist_k .* h);
    
    dist_next = zeros(Nk, 1);
    
    for ik= 1:Nk
        
        mass_here= dist_k(ik);
        
        replacing= mass_here * h(ik);
        
        staying= mass_here * (1 - h(ik));
        
        dist_next(1)= dist_next(1) + replacing;
        
        if ik < Nk
            ik_next= ik + 1;
        else
            ik_next= Nk;
        end
        dist_next(ik_next)= dist_next(ik_next) + staying;
        
    end
    
    dist_next= dist_next / sum(dist_next);
    
    dist_k= dist_next;
end

figure;
plot(1:T_agg, I_path, '-o');
xlabel('Time');
ylabel('Aggregate replacement rate');
title('Aggregate investment with fixed A (Figure 3-style)');

%%Figure 3 results:

% With A fixed, the model already produces lumpy investment at the aggregate level.
% Starting from a uniform distribution of capital ages, the aggregate replacement rate
% adjusts over time as the cross-section of firms moves toward its stationary
% distribution. Even though the hazard of replacement is purely age-dependent (given A),
% the interaction between the age profile of capital and the replacement rule generates
% spikes in aggregate investment when many firms hit the replacement region around the
% same time.
%
% In the long run, the aggregate replacement rate settles down to a stationary level:
% the flow of firms replacing each period stabilizes because the cross-sectional
% distribution of capital ages has converged. This shows that non-convex adjustment
% costs and (S,s)-type behavior can create aggregate lumpy investment even when the
% aggregate state A is constant and there are no aggregate shocks.



T_agg_MC = 200;

I_path_MC= zeros(T_agg_MC,1);
A_index_MC= zeros(T_agg_MC,1);
A_series_MC= zeros(T_agg_MC,1);

 



dist_k= ones(Nk, 1) / Nk;

A_index_MC(1)= 1;
A_series_MC(1)= A_grid(A_index_MC(1));

for t= 1:T_agg_MC
    
    ia= A_index_MC(t);
    h_t = hazard(:,ia);
    
   
    I_path_MC(t)= sum(dist_k .* h_t);
    
    dist_next= zeros(Nk, 1);
    
    for ik= 1:Nk
        
        mass_here= dist_k(ik);
        
        replacing= mass_here * h_t(ik);
        
        staying=mass_here * (1 - h_t(ik));
        
        dist_next(1)= dist_next(1) + replacing;
        
        if ik < Nk
            ik_next= ik + 1;
        else
            ik_next= Nk;
        end
        dist_next(ik_next)= dist_next(ik_next) + staying;
        
    end
    
    
    dist_next= dist_next / sum(dist_next);
    dist_k= dist_next;
    
    if t < T_agg_MC
        u= rand;
        curA= A_index_MC(t);
        
        if u <= P_A(curA, 1)
            A_index_MC(t+1)= 1;
        else
            A_index_MC(t+1)= 2;
        end
        
        A_series_MC(t+1)= A_grid(A_index_MC(t+1));
    end
    
end

figure;

subplot(2,1,1);
plot(1:T_agg_MC, I_path_MC, '-');
xlabel('Time');
ylabel('Aggregate replacement rate');
title('Aggregate investment with Markov A (Figure 4-style)');

subplot(2,1,2);
stairs(1:T_agg_MC, A_series_MC, '-');
xlabel('Time');
ylabel('Aggregate productivity A_t');
title('Aggregate state over time');
%% Figure results:

% When A follows the Markov process, the hazard of replacement becomes state-dependent:
% for a given capital age, the probability of replacement depends on whether A is high
% or low. In the simulation, this makes the aggregate replacement rate comove with the
% aggregate state: investment tends to be different in high-A vs low-A periods.
%
% The interaction of (i) the time-varying hazard (through A) and (ii) the evolving
% cross-sectional distribution of capital ages generates fluctuations in aggregate
% investment that are both uneven and history-dependent. Depending on the current
% distribution of ages, a shift in A can trigger a wave of replacements or have only
% a modest effect. This matches the idea in Cooper, Haltiwanger, and Power (1999)
% that aggregate investment can show large "lumps and bumps" even when the underlying
% shocks are relatively simple.
%
% Overall, the model implies that firm investment behavior is highly non-smooth at
% both the micro and macro level. At the firm level, machines are kept in place for
% several periods and then replaced in a discrete jump when the state hits the
% replacement region. At the aggregate level, the timing and size of investment spikes
% depend on the joint distribution of capital ages and the history of aggregate shocks.
