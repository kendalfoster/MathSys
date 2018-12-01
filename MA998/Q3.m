clear variables
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',14)

Nb = 1024;
b_vec = linspace(-1,1,Nb);
a_vec = linspace(1,8,8);

Ntrans = 256;
NT = 256;
xstart = sqrt(2)/10;
lambda = zeros(1,Nb);

xtrans = zeros(Ntrans,Nb);
x = zeros(NT,Nb);

%for i = 1:8
%a = a_vec(i);
a= 6;

for k = 1:Nb   
    xtrans(1,k) = xstart;
    b = b_vec(k); 
    
    
    for j = 1:Ntrans-1
        xtrans(j+1,k) = exp(-a*(xtrans(j,k)^2))+b;
    end    
    
    x(1,k) = xtrans(Ntrans,k);    
    
    for q = 1:NT-1
        x(q+1,k)= exp(-a*(x(q,k)^2))+b;
    end
    
    
    lambda(1,k) = (sum( log(-2.*a.*xtrans(:,k).*exp(-a.*(xtrans(:,k).^2)) ))+sum( log(-2.*a.*x(:,k).*exp(-a.*(x(:,k).^2)) )) )*(1/(NT+Ntrans));

end

zero = zeros(Nb);
 
figure(1)
subplot(2,1,1)
plot(b_vec,x,'k.','MarkerSize',1)
xlabel('Control parameter b')
ylabel('lim_{n\rightarrow\infty} x_n')
title ('Bifurcation Diagram')

subplot(2,1,2)
plot(b_vec,lambda,'r',b_vec,zero,'b')
xlabel('Control parameter b')
ylabel('\lambda')
title ('Lyapunov Exponent')

%end