no_of_independent_trials = 20;

theta = 1.01;

for itr=1:no_of_independent_trials
    
    clc;
    disp(['Independent Trial No: ',num2str(itr)])
    
    no_of_input = 500000;
    
    input=rand(1,no_of_input) - 0.5;
    
    N=10;
        
    AEFLN_order = 1;
    
    input_buffer= zeros(1,N);
        
    M = (2*AEFLN_order+1)*N + 1;
    
    AEFLN_weights= zeros(1,M);
    
    mu_weight= 0.001;
    
    mu_a= 0.05;
        
    noise = awgn(input,30)-input;
    
    a(1)= 0;
    SNR = 30;

    g1 = zeros(1,no_of_input); 
    g2 = zeros(1,no_of_input); 
    g3 = zeros(1,no_of_input);

    for i=3:no_of_input
        g1(i) = exp(0.5*input(i))*(sin(pi*input(i)) + 0.3*sin(3*pi*input(i-2)) + 0.1*sin(5*pi*input(i)));
    end
    g1 = awgn(g1, SNR);

    for i=1:no_of_input
        q = (3/2)*input(i) - (3/10)*input(i)^2;
        if q > 0, rho=4; else rho=0.5; end
        g2(i) = 2*((1/(1+exp(-rho*q))) - 0.5);
    end
    g2 = awgn(g2, SNR);

    chi = 0.1;
    for i=1:no_of_input
       if abs(input(i))>=0 && abs(input(i))<chi
           g3(i) = (2/(3*chi))*input(i);
       elseif abs(input(i))>=chi && abs(input(i))<(2*chi)
           g3(i) = sign(input(i))*(3 - (2 - abs(input(i))/chi)^2)/3;
       else
           g3(i) = sign(input(i));   
       end
    end
    g3 = awgn(g3, SNR);

    g = zeros(1,no_of_input);
    
    for i=1:length(input)
        
        emp = [];
        for mn = 1:N
            emp = [emp, theta^(N-mn)];
        end

        input_buffer=[input(i) input_buffer(1:end-1)];
        input_buffer = (diag(emp)*input_buffer')';
%         theta_mat = diag([1, 1, emp, 1, emp, 1, emp]);
        
        AE_FEB=[];
        for k =1:N
            for l =1:AEFLN_order
                AE_FEB=[AE_FEB, exp(-1*a(i)*abs(input_buffer(k)))*sin(pi*l*input_buffer(k)), exp(-1*a(i)*abs(input_buffer(k)))*cos(pi*l*input_buffer(k))];
            end
        end
        
        AEFLN_after_FEB= ([1,input_buffer,AE_FEB]')';
        
        if i < no_of_input/3
            g(i)=g2(i);
        elseif i >= no_of_input/3 && i < 2*no_of_input/3
            g(i)=g3(i);
        else
            g(i)=g1(i);
        end
        desired_output(i) = g(i);
    
        AEFLN_output(i)= AEFLN_weights * AEFLN_after_FEB';
    
        error(i)= desired_output(i) - AEFLN_output(i);
    
    
        z=[];
        for k =1:N
            for l =1:AEFLN_order
                z=[z, -1*abs(input_buffer(k))*exp(-1*a(i)*abs(input_buffer(k)))*sin(pi*l*input_buffer(k)), -1*abs(input_buffer(k))*exp(-1*a(i)*abs(input_buffer(k)))*cos(pi*l*input_buffer(k))];
            end
        end
        
        z_final= ([0,zeros(1,N),z]')';   
        
    
        a(i+1) = a(i) + mu_a * error(i) * z_final * AEFLN_weights';
            
        AEFLN_weights= AEFLN_weights + mu_weight * error(i) * AEFLN_after_FEB;
        
    end
    err(itr,:)=error.^2;
end
%%
length_of_smoothing_filter = 1000;

smoothing_filter_coeff = (1/length_of_smoothing_filter)*ones(1,length_of_smoothing_filter);
for i=1:no_of_independent_trials
    err_smooth(i,:) = filter(smoothing_filter_coeff,1,err(i,:));
end
% figure;
plot(10*log10(mean(err_smooth)));
