no_of_input = 50000;
N = 10;                   
fln_order = 2;               
M = (2*fln_order+1)*N + 1;

lambda = 0.999;        
delta = 1e-3;         
mu_a = 0.0005;              
theta = 1.003;

input = rand(1,no_of_input) - 0.5;
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

Wb = zeros(M,1);         
P = (1/delta) * eye(M);  

a(1) = 0;                         
G = 0;                           

x_buffer = zeros(1,N);
yhat = zeros(1,no_of_input);
err = zeros(1,no_of_input);

for i = 1:no_of_input
    
    x_buffer = [input(i), x_buffer(1:end-1)];
        
    fln_input = [1, x_buffer, exp(-a(i) * abs(x_buffer)) .* sin(pi * x_buffer), exp(-a(i) * abs(x_buffer)) .* cos(pi * x_buffer), exp(-a(i) * abs(x_buffer)) .* sin(2*pi * x_buffer), exp(-a(i) * abs(x_buffer)) .* cos(2*pi * x_buffer)]';   
    if i < no_of_input/3
        g(i)=g2(i);
    elseif i >= no_of_input/3 && i < 2*no_of_input/3
        g(i)=g3(i);
    else
        g(i)=g1(i);
    end
    yhat(i) = Wb' * fln_input;
    err(i) = g(i) - yhat(i);
    
    Tn = theta^(i-1);
    Pz = P * fln_input;                    
    denom = lambda + Tn^2*(fln_input' * Pz);  
    K = Pz / denom;               
    Wb = Wb + Tn^2*K * err(i);          
    P = (1 / lambda) * (P - Tn^2*(Pz * Pz') / denom);  
    
    z_final = [0, zeros(1,N), (-abs(x_buffer) .* exp(-a(i) * abs(x_buffer))) .* sin(pi * x_buffer), (-abs(x_buffer) .* exp(-a(i) * abs(x_buffer))) .* cos(pi * x_buffer), (-abs(x_buffer) .* exp(-a(i) * abs(x_buffer))) .* sin(2*pi * x_buffer), (-abs(x_buffer) .* exp(-a(i) * abs(x_buffer))) .* cos(2*pi * x_buffer)]';  
    
    
    h_n = err(i) * (Wb' * z_final);   
    
    G = lambda * G + h_n;
    
    a(i+1) = a(i) + mu_a * G;
    
end

N = 200;
aef_err = err.^2;
smooth_err = filter(ones(1,N)/N,1,aef_err);

% figure;
plot(10*log10(smooth_err)); 
