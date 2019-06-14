clc; close all; pause(1);
clear; 
f1 = 10;
f2 = 12;
f3 = 11;
Vout = 4;              % 4 for sensor observations & 7 for compressed level
Const = 'v4';          
L0 = 35;               
Lend = 40;             
L = 5000;             
Nlevels0 = 3;          
n_std = 0.3;           
                       
self = 1;                % Using previous step as reference in calculation of error
nFlag = 1;             % 0 for norm-2 and 1 for norm-1
LFlag = 1;             % 0 for just the immediate level-based construction & o for all so-far selected marginal 

delta = 0.2;          % delta changed from 0.1 
MA =10;               % moving average window size
perc = 5/100;
step = 0.01;
bins = 500;
Thr = 0.2;
%Mu = 1;
mu = 1;
dyn_end = 200;
Navg = 5;             % maximum 50
jmp = 1;
Npnt = fix(32/jmp);

Error_opt = zeros(1, Npnt);
Error_prop = zeros(1, Npnt);
Error_unif = zeros(1, Npnt);
Error_learn_prop = zeros(1, Npnt);
Error_learn_unif = zeros(1, Npnt);
Delta_prop_avg = zeros(1, Npnt);
Delta_unif_avg = zeros(1, Npnt);
Slevels = zeros(1, Npnt);
Sensors = zeros(L,7+MA+2);

Cost_opt = zeros(size(Slevels));
Cost_prop = zeros(size(Slevels));
Cost_unif = zeros(size(Slevels));
[XI,YI] = meshgrid(0:1:100);

Peak = 100;
Sensors(:,1:2) = Peak * rand(L,2);
tic;
% ======================================================= Proposed
for rept = 1 : Navg
    error_1 = 0;
    error_2 = 0;
    W = zeros(1,5);
    L_bank = zeros(Nlevels0 + Npnt + 2, Npnt);
    Nlevels = Nlevels0;
    Delta = delta;
    Levels = linspace(L0,Lend,Nlevels);
    D_lock = 0;                    
    
    % Generation of the large scale Gaussians
    name = sprintf('NewWave_%d',rept);
    load(name);
    %Sensors(:,1:2) = Peak * rand(L,2);
    X = Sensors(:,1);
    Y = Sensors(:,2);
    Z = zeros(L,1);
    for k = 1: N
        Z = Z + 2*Peak * (P(k,3) - 0.5) / sqrt(2*pi) / Sigma * exp(-((X-P(k,1)).^2 + (Y-P(k,2)).^2)/2/Sigma^2);
    end,

    % Genmeration of small scale Gaussians
    for k = 1 : Ns
        Z = Z + 5 * (Ps(k,3) - 0.5) / sqrt(2*pi) / Ss * exp(-((X-Ps(k,1)).^2 + (Y-Ps(k,2)).^2)/2/Ss^2);
    end,

    % Snapshot normalization
    Z = Peak / 2 * Z / (max(max(Z)) - min(min(Z)))+50;
    Sensors(:,3) = Z;        
    Sensors(:,4) = Z + n_std * randn(size(Z));   
    
    Sensors(:,5:7) = 0; 
    for pnt = 1 : Npnt   
%{        
        if (pnt < 6)
            jmp = 2;
        elseif (pnt >= 6 && pnt < 15)
            jpm = 3;
        elseif (pnt >= 15)
            jmp = 5;
        end,
%}
        Sensors(:,6:7) = 0;         
        for clevel = 1 : length(Levels)
            index = find(abs(Sensors(:,4) - Levels(clevel)) <= Delta);        
            Sensors(index, 6) = 1;
            Sensors(index, 7) = Levels(clevel);            
        end,
        index = find(Sensors(:,6) == 1);
%{        
        if (LFlag == 1)
            Filtered = sum(abs(Sensors(index, 6) - Sensors(index, 5)));
        elseif (LFlag == 0)
            Filtered =  length(index);
        end,
%}        
        Filtered = sum(abs(Sensors(index, 6) - Sensors(index, 5)));
        Sensors(index, 5) = 1;
        
        Cost_prop(pnt) = Cost_prop(pnt) + Filtered;

        if (LFlag == 1)        
            Index_select = find(Sensors(:,5) == 1);
        elseif (LFlag == 0)
            Index_select = find(Sensors(:,6) == 1);
        end
        
        X0 = Sensors(Index_select, 1);
        Y0 = Sensors(Index_select, 2);
        Z0 = Sensors(Index_select, Vout);

        ZI = griddata(X0,Y0,Z0,XI,YI,Const);          % linear   nearest      cubic     v4(biharmonic interpolation)
        if (nFlag == 0)
            error_1 = sqrt(sum(sum((ZI - Data).^2)/(101*101)));
        elseif (nFlag == 1)
            error_1 = sum(sum(abs(ZI - Data))/(101*101));
        end
        Error_prop(pnt) = Error_prop(pnt) + error_1;
        
        if (pnt > 1)
           %Corr =  corr2(Data_old,ZI);
           %Error_learn_prop(pnt) = Corr;
           if (nFlag == 0)
                error_2 = sqrt(sum(sum((ZI - Data_old).^2)/(101*101)));
           elseif (nFlag == 1)
                error_2 = sum(sum(abs(ZI - Data_old))/(101*101));
           end
           Error_learn_prop(pnt) = Error_learn_prop(pnt) + error_2;       
        end,
        Data_old = ZI; 
        
        L_bank(1,pnt) = pnt+Nlevels0-1;
        L_bank(2,pnt) = Delta;
        L_bank(3:2+Nlevels, pnt) = Levels;

        Slevel(pnt) = Nlevels;    
        [p,xx] = ksdensity(ZI(:),'npoints',bins);
        %[p,xx] = hist(ZI(:),bins);
        %p = p / sum(p);
        Nlevels = Nlevels + jmp;
        Levels = Lloyd_Max_2(p, xx, Nlevels);        
        
        if (self == 0)
            error = error_1;
            if (pnt ==1)
                error_old = error;
            end,
        end,        
        if (self == 1)
            error = error_2;
            if (pnt <= 2)
                error_old = error;
            end,
        end,
        %mu = Mu * (Npnt-pnt)/Npnt;
        %mu = Mu * (pnt)/Npnt;
        if (pnt > 1)
            Delta = Delta *(1 + mu*(error - error_old)/(error + error_old));
            %Delta = Delta *(1+mu*(Error_prop(pnt) - Error_prop(pnt-1))/(Error_prop(pnt) + Error_prop(pnt-1)));            
        end,
        error_old = error;
    end,    
    Delta_prop_avg = Delta_prop_avg + L_bank(2,:);
end,
Delta_done = Delta;
Cost_prop = Cost_prop / Navg;
Error_prop = Error_prop / Navg;
Delta_prop_avg = Delta_prop_avg / Navg;
Error_learn_prop = Error_learn_prop / Navg;
% figure(2); mesh(XI,YI,ZI); pause(1);

%{
[Min, indx] = min(Error_prop);
ix = L_bank(1,indx);
Levels = L_bank(3:2+ix , indx);
Levels_s = Levels;
Delta = L_bank(2,indx);
%}
% ===========================================
% Dynamic case
% ===========================================
x_s = step;     % x direction shift increment
y_s = step;     % y direction shift increment
reps = 0;
Sensors(:,end) = Sensors(:,4);
Cost = zeros(ceil(dyn_end/MA),1);
Error = zeros(ceil(dyn_end/MA),1);
m = 1;
Error_0 = Error_prop(end);
xx0 = 0;
yy0 = 0;
Delta = Delta_done
for dyn = 1:dyn_end
% Generation of the large scale Gaussians
    X = Sensors(:,1);
    Y = Sensors(:,2);
    %Sensors(:,6:7) = 0;
    Z = zeros(L,1);
    Data_new = zeros(size(XI));
    
    if (abs(cos(2*pi*dyn*step)) > Thr)
        xx0 = xx0 + x_s;
        yy0 = yy0 + y_s;
        reps = reps + 1;
    end,
    
    for k = 1: N
        Z = Z + 2*Peak * (P(k,3) - 0.5) / sqrt(2*pi) / Sigma * exp(-((X-P(k,1)-xx0).^2 + (Y-P(k,2)-yy0).^2)/2/Sigma^2);
        Data_new = Data_new + 2*Peak * (P(k,3) - 0.5) / sqrt(2*pi) / Sigma * exp(-((XI-P(k,1)-xx0).^2 + (YI-P(k,2)-yy0).^2)/2/Sigma^2);
    end,
    
% Generation of small scale Gaussians
    for k = 1 : Ns
        Z = Z + 5 * (Ps(k,3) - 0.5) / sqrt(2*pi) / Ss * exp(-((X-Ps(k,1)).^2 + (Y-Ps(k,2)).^2)/2/Ss^2);
        Data_new = Data_new + 5 * (Ps(k,3) - 0.5) / sqrt(2*pi) / Ss * exp(-((XI-Ps(k,1)).^2 + (YI-Ps(k,2)).^2)/2/Ss^2);
    end,
    Z = Peak / 2 * Z / (max(Z) - min(Z))+50;
    Data_new = Peak / 2 * Data_new / (max(max(Data_new)) - min(min(Data_new)))+50;
    
    Sensors(:,4) = Z + n_std * randn(size(Z));     
    %Sensors(:,7+1:7+MA-1) = Sensors(:,7+2:7+MA);
    %Sensors(:,7+MA) = Sensors(:,4);
    Sensors(:,5:7) = 0;
    
    if (dyn > MA)
            %Sensors(:,end-1) = mean(Sensors(:,7+1:7+MA)');
            %diss = find(abs(Sensors(:,end) - Sensors(:,end-1)) > 2*Delta);        
        %if (length(diss)/L > perc) && (rem(dyn,MA) == 0)
        if (rem(dyn,MA) == 0)
            %Sensors(:, 5:7) = 0;
            for clevel = 1 : length(Levels)
                index = find(abs(Sensors(:,4) - Levels(clevel)) <= Delta);  
                Sensors(index, 6) = 1;
                %Sensors(index, 7) = Levels(clevel);                    
            end,
            Index = find(Sensors(:, 6) == 1);
            Cost(m) = length(Index); %+ length(diss);

            X1 = Sensors(Index, 1);
            Y1 = Sensors(Index, 2);
            Z1 = Sensors(Index, 4);
            Zx = griddata(X1,Y1,Z1,XI,YI,Const);          % linear   nearest    cubic     v4(biharmonic interpolation)
            if (nFlag == 0)
                Error(m) = sqrt(sum(sum((Zx - Data_new).^2)/(101*101)));  
            elseif (nFlag == 1)
                Error(m) = sum(sum(abs(Zx - Data_new))/(101*101));
            end,
            [p,xx] = ksdensity(Zx(:),'npoints',bins);
            Levels = Lloyd_Max_2(p, xx, length(Levels));  
            %Error_0 = Error(m);
            %Sensors(:,end) = Sensors(:,end-1);            
        end,
%{
        if (length(diss)/L <= perc)
            Cost(m) = length(diss);
            Error(m) = Error_0;
            reps = reps + 1;
        end, 
%}        
        m = m+1;                
    end,    
end,

%figure(5); mesh(XI,YI,Data_new-Zx);title('Error Magnitude');
%figure(6); mesh(XI,YI,Data_new);
%figure(5); hold on; plot3(Sensors(:,1),Sensors(:,2),Z,'.');
%figure(6); hold on; mesh(XI,YI,Zx);
figure(101); plot(20*log10(Error),'r*'); grid;title('performance'); ylabel('MSE in dB'); xlabel('update instant');
figure(100); plot(Cost/L*100,'bs'); grid;title('cost'); ylabel('Percentage of sensors'); xlabel('update instant');
toc/60
%reps
pause(1);

tic;
% ================================================================= Uniform Adapted
%Sensors(:,1:2) = Peak * rand(L,2);
for rept = 1 : Navg
    error_1 = 0;
    error_2 = 0;
    W = zeros(1,5);
    LU_bank = zeros(Nlevels0 + Npnt + 2, Npnt);    
    Nlevels = Nlevels0;
    Delta = delta;
    %Levels = linspace(L0,Lend,Nlevels);
    Levels = linspace(L0,Lend,Nlevels+2);
    Levels = Levels(2:end-1);
    D_lock = 0;     
    
    % Generation of the large scale Gaussians
    name = sprintf('NewWave_%d',rept);
    load(name);
    %Sensors(:,1:2) = Peak * rand(L,2);
    X = Sensors(:,1);
    Y = Sensors(:,2);
    Z = zeros(L,1);
    for k = 1: N
        Z = Z + 2*Peak * (P(k,3) - 0.5) / sqrt(2*pi) / Sigma * exp(-((X-P(k,1)).^2 + (Y-P(k,2)).^2)/2/Sigma^2);
    end,

    % Genmeration of small scale Gaussians
    for k = 1 : Ns
        Z = Z + 5 * (Ps(k,3) - 0.5) / sqrt(2*pi) / Ss * exp(-((X-Ps(k,1)).^2 + (Y-Ps(k,2)).^2)/2/Ss^2);
    end,

    % Snapshot normalization
    Z = Peak / 2 * Z / (max(max(Z)) - min(min(Z)))+50;
    Sensors(:,3) = Z;        
    Sensors(:,4) = Z + n_std * randn(size(Z));   
    
    Sensors(:,5:7) = 0; 
    for pnt = 1 : Npnt
%{
        if (pnt < 6)
            jmp = 2;
        elseif (pnt >= 6 && pnt < 15)
            jpm = 3;
        elseif (pnt >= 15)
            jmp = 5;
        end,
%}        
        Sensors(:,6:7) = 0;         
        for clevel = 1 : length(Levels)
            index = find(abs(Sensors(:,4) - Levels(clevel)) <= Delta);        
            Sensors(index, 6) = 1;
            Sensors(index, 7) = Levels(clevel);            
        end,
        index = find(Sensors(:,6) == 1);        
%{        
        if (LFlag == 1)
            Filtered = sum(abs(Sensors(index, 6) - Sensors(index, 5)));
        elseif (LFlag == 0)
            Filtered =  length(index);
        end,
%}        
        Filtered = sum(abs(Sensors(index, 6) - Sensors(index, 5)));
        Sensors(index, 5) = 1;
        
        Cost_unif(pnt) = Cost_unif(pnt) + Filtered;

        if (LFlag == 1)        
            Index_select = find(Sensors(:,5) == 1);
        elseif (LFlag == 0)
            Index_select = find(Sensors(:,6) == 1);
        end
        X0 = Sensors(Index_select, 1);
        Y0 = Sensors(Index_select, 2);
        Z0 = Sensors(Index_select, Vout);

        ZI = griddata(X0,Y0,Z0,XI,YI,Const);          % linear   nearest      cubic     v4(biharmonic interpolation)
        if (nFlag == 0) 
            error_1 = sqrt(sum(sum((ZI - Data).^2)/(101*101)));
        elseif (nFlag == 1)
            error_1 = sum(sum(abs(ZI - Data))/(101*101));
        end,
        Error_unif(pnt) = Error_unif(pnt) + error_1;
        
        if (pnt > 1)
           %Corr =  corr2(Data_old,ZI);
           %Error_learn_unif(pnt) = Corr;
           if (nFlag == 0)
                error_2 = sqrt(sum(sum((ZI - Data_old).^2)/(101*101)));
           elseif (nFlag == 1)
                error_2 = sum(sum(abs(ZI - Data_old))/(101*101));
           end
           Error_learn_unif(pnt) = Error_learn_unif(pnt) + error_2;       
        end,
        Data_old = ZI; 
        
        LU_bank(1,pnt) = pnt+Nlevels0-1;
        LU_bank(2,pnt) = Delta;
        LU_bank(3:2+Nlevels, pnt) = Levels;

        Nlevels = Nlevels + jmp;
        L_0 = min(min(ZI));
        L_end = max(max(ZI));
        %Levels = linspace(L_0,L_end,Nlevels);
        Levels = linspace(L_0,L_end,Nlevels+2);
        Levels = Levels(2:end-1);

        if (self == 0)
            error = error_1;
            if (pnt ==1)
                error_old = error;
            end,
        end,        
        if (self == 1)
            error = error_2;
            if (pnt <= 2)
                error_old = error;
            end,
        end,
        %mu = Mu * (Npnt - pnt)/Npnt;
        %mu = Mu * (pnt)/Npnt;
        if (pnt > 1)
            Delta = Delta *(1 + mu*(error - error_old)/(error + error_old));
            %Delta = Delta *(1+mu*(Error_unif(pnt) - Error_unif(pnt-1))/(Error_unif(pnt) + Error_unif(pnt-1)));            
        end,
        error_old = error;
    end,
    Delta_unif_avg = Delta_unif_avg + LU_bank(2,:);
end,    

Cost_unif = Cost_unif / Navg;
Error_unif = Error_unif / Navg;
Delta_unif_avg = Delta_unif_avg / Navg;
Error_learn_unif = Error_learn_unif / Navg;
toc/60

figure(21); plot(Slevel, L_bank(2,:),'r', Slevel, LU_bank(2,:),'g');grid;
figure(22); plot(Slevel, Delta_prop_avg,'r', Slevel, Delta_unif_avg,'g');grid;

Error_learn_prop(1) = Error_learn_prop(2);
Error_learn_unif(1) = Error_learn_unif(2);
figure(23); plot(Slevel, 20*log10(Error_learn_prop),'-dr', Slevel, 20*log10(Error_learn_unif),'-sg');grid on;
pause(1);

tic;
% ======================================================================== Optimal
%Sensors(:,1:2) = Peak * rand(L,2);
for rept = 1 : Navg
    %[p,xx] = hist(Data(:),bins);
    %p = p / sum(p);
    [p,xx] = ksdensity(Data(:),'npoints',bins);
    Nlevels = Nlevels0;
    Delta = delta;
    D_lock = 0;        
    
    % Generation of the large scale Gaussians
    name = sprintf('NewWave_%d',rept);
    load(name);
    %Sensors(:,1:2) = Peak * rand(L,2);
    X = Sensors(:,1);
    Y = Sensors(:,2);
    Z = zeros(L,1);
    for k = 1: N
        Z = Z + 2*Peak * (P(k,3) - 0.5) / sqrt(2*pi) / Sigma * exp(-((X-P(k,1)).^2 + (Y-P(k,2)).^2)/2/Sigma^2);
    end,

    % Genmeration of small scale Gaussians
    for k = 1 : Ns
        Z = Z + 5 * (Ps(k,3) - 0.5) / sqrt(2*pi) / Ss * exp(-((X-Ps(k,1)).^2 + (Y-Ps(k,2)).^2)/2/Ss^2);
    end,

    % Snapshot normalization
    Z = Peak / 2 * Z / (max(max(Z)) - min(min(Z)))+50;
    Sensors(:,3) = Z;        
    Sensors(:,4) = Z + n_std * randn(size(Z));   
    Sensors(:,5:7) = 0;
    for pnt = 1 : Npnt
%{        
        if (pnt < 6)
            jmp = 2;
        elseif (pnt >= 6 && pnt < 15)
            jpm = 3;
        elseif (pnt >= 15)
            jmp = 5;
        end,     
%}        
        Levels = Lloyd_Max_2(p, xx, Nlevels);
        Sensors(:,6:7) = 0;         
        for clevel = 1 : length(Levels)
            index = find(abs(Sensors(:,4) - Levels(clevel)) <= Delta);        
            Sensors(index, 6) = 1;
            Sensors(index, 7) = Levels(clevel);            
        end,
        index = find(Sensors(:,6) == 1);        
%{        
        if (LFlag == 1)
            Filtered = sum(abs(Sensors(index, 6) - Sensors(index, 5)));
        elseif (LFlag == 0)
            Filtered =  length(index);
        end,
%}        
        Filtered = sum(abs(Sensors(index, 6) - Sensors(index, 5)));
        Sensors(index, 5) = 1;
        
        Cost_opt(pnt) = Cost_opt(pnt) + Filtered;

        if (LFlag == 1)        
            Index_select = find(Sensors(:,5) == 1);
        elseif (LFlag == 0)
            Index_select = find(Sensors(:,6) == 1);
        end
        X0 = Sensors(Index_select, 1);
        Y0 = Sensors(Index_select, 2);
        Z0 = Sensors(Index_select, Vout);

        ZI = griddata(X0,Y0,Z0,XI,YI,Const);          % linear   nearest      cubic     v4(biharmonic interpolation)
        if (nFlag == 0)
            error = sqrt(sum(sum((ZI - Data).^2)/(101*101)));
        elseif (nFlag == 1)
            error = sum(sum(abs(ZI - Data))/(101*101));
        end
        Error_opt(pnt) = Error_opt(pnt) + error;

        Nlevels = Nlevels + jmp;
    end,
end,

Cost_opt = Cost_opt / Navg;
Error_opt = Error_opt / Navg;

toc/60
Cost = [Cost_prop ; Cost_unif ; Cost_opt];
for kk = 2 : length(Slevel)
    Cost(: , kk) = Cost(: , kk) + Cost(: , kk-1);
end,
Cost = Cost/L*100;

% ========================================================================
figure(f1); plot(Slevel(1:end), 20*log10(Error_prop(1:end)),'-sr', Slevel(1:end), 20*log10(Error_unif(1:end)),'-dg', Slevel(1:end), 20*log10(Error_opt(1:end)), '-ob'); grid on;ylabel('Mean Error in dB'); xlabel('# of levels');
figure(f3); plot(Slevel(1:end), (Error_prop(1:end)),'-sr', Slevel(1:end), (Error_unif(1:end)),'-dg', Slevel(1:end), (Error_opt(1:end)), '-ob'); grid on;ylabel('Mean Error'); xlabel('# of levels');
figure(f2); plot(Slevel(1:end), Cost_prop(1:end)/L*100,'r', Slevel(1:end), Cost_unif(1:end)/L*100,'g', Slevel(1:end), Cost_opt(1:end)/L*100,'b'); grid on; ylabel('Percentage of sensors'); xlabel('# of levels');
%figure(21); plot(Slevel, L_bank(2,:),'r', Slevel, LU_bank(2,:),'g');grid;
%figure(22); plot(Slevel, Delta_prop_avg,'r', Slevel, Delta_unif_avg,'g');grid;
figure(24); plot(Slevel, Cost(1,:),'r', Slevel, Cost(2,:),'g', Slevel, Cost(3,:),'b'); grid on; ylabel('Percentage of sensors'); xlabel('# of levels');

%Error_learn_prop(1) = Error_learn_prop(2);
%Error_learn_unif(1) = Error_learn_unif(2);
%figure(23); plot(Slevel, 20*log10(Error_learn_prop),'-dr', Slevel, 20*log10(Error_learn_unif),'-sg');grid on;

Relative_costs = [sum(Cost_opt)/L   sum(Cost_prop)/L    sum(Cost_unif)/L]
pause(1);

%{
Z_all = griddata(Sensors(:,1),Sensors(:,2),Sensors(:,4),XI,YI,'v4');          % linear   nearest      cubic     v4(biharmonic interpolation)
if (nFlag == 0)
    error_v4_all = sqrt(sum(sum((Z_all - Data).^2)/(101*101)));
elseif (nFlag == 1)
    error_v4_all = sum(sum(abs(Z_all - Data))/(101*101));
end
All_Sensor_Error = [20*log10(error_v4_all)]

figure(1); plot(Slevel(2:end), 20*log10(Error_learn_prop(1:end-1)),'r', Slevel(2:end), 20*log10(Error_learn_unif(1:end-1)),'g'); grid on;
%}