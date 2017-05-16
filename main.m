% Cheikh Toure

clear all

lw=2; %Linewidth
fs=18; %Fontsize
fw='Bold'; %FontWeight
fsa=16; %Fontsize

tab_n = 5:1:15;
n=length(tab_n);

tab_s = 1 + tab_n; 

Tn = zeros(size(tab_n));
Ts = zeros(size(tab_n));


for i = 1:n
    
    N = 2^(tab_n(i));
    s = tab_s(i);
    
    x = rand(N,1);
    
   idx = randi(N,s,1);
    
    tic
    hadamardn(x);
    toc
    Tn(i) = toc
    
    tic
    hadamards(x,idx);
    toc
    Ts(i) = toc
end

 plot(tab_n,Tn,tab_n,Ts,'Linewidth',lw);
 xlabel('Absciss (in log_2(n) where n is the length of the input)','FontSize',fs,'FontWeight',fw);
 ylabel('Ordinate (time in sec)','FontSize',fs,'FontWeight',fw);
 title('Time complexities','FontSize',fs,'FontWeight',fw);
 set(gca,'FontWeight',fw,'FontSize',fsa);
 legend({'Time (hadamardn)','Time (hadamard with s = log_2(n))'},'FontWeight',fw,'FontSize',fsa);
