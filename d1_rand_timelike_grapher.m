% constants
tic;
L = 100;
fill = [0.5];

% basic time step (brillouin zone lattice spacing) 
u = 1/4;
tau0 = 2*pi/(4*u);

% fixed tau, Krange number of points, K indexes up to Krange
tau = 2*tau0;
Kmax = 60;
K=10:Kmax;

% calculate the normalizations
uv = 1e-16;
dos = @(E) 1./(E.*log(1./E).^3);

totaldensity = @(c) integral(@(E) c.*dos(E), uv, 1/2, 'ArrayValued',true);
norm = fzero(@(x) totaldensity(x)-1/2, [0,2]);

% store the data: (fill, length(K), rep rate)
data = zeros(length(fill), length(K));

% fetch the random diagonalized hamiltonian from the helper method
for i = 1:length(fill)
    for j = 1:length(K)
        data(i,j) = d1_rand_timelike_dos(L, K(j), tau, uv, norm);
    end
end

S_mean = mean(data, 3);
toc

% slope = zeros(2,length(fill));
% kf = zeros(1, length(fill));
% for i = 1:length(fill)
%     X = [ones(length(K),1) log(log(K)).']; 
%     slope(:,i) = X\(S_mean(i,:).');
% end
% disp(slope(2,:))

Kmin = 1;

c = zeros(2,3);

modelfun = @(b,x)(b(1)*log(log(x)-log(b(2))));

beta0 = [0.5;1];
[beta,R,J,CovB,MSE] = nlinfit(((Kmin+9):(Kmax))*tau,S_mean(1,:),modelfun,beta0);
beta = real(beta);
c(1,:) = [beta(1) beta(2) MSE];

modelfun = @(b,x)(b(1)*(log(x)-log(b(2))));

beta0 = [0.5;1];
[beta2,R2,J2,CovB2,MSE2] = nlinfit(((Kmin+9):(Kmax))*tau,S_mean(1,:),modelfun,beta0);
beta2 = real(beta2);
c(2,:) = [beta2(1) beta2(2) MSE2];

disp(c)

colorstring = 'kbgry';
close all;
figure; cla;
hold on
for i = 1:length(fill)
    % (1/3)*log((L/pi).*sin(pi.*lx./L))
    % scatter((1/3)*log((L/pi).*sin(pi.*lx./L)),S_mean(i,:),8,'Color',colorstring(i));
    scatter(log(log(K*tau)),S_mean(i,:),8,'Color',colorstring(i));
    %scatter(log(K*tau),S_mean(i,:),8,'Color',colorstring(i));
end
for i=1:length(fill)
    plot(log(log(K*tau)),X*slope(:,i),'Color',colorstring(i));
    %plot(log(K*tau),X*slope(:,i),'Color',colorstring(i));
end

title('L= 100, t = [10,100], tau = 2tau0, T=0')
xlabel('log(log(t))') 
ylabel('S_A') 
%lgd = legend('n=0.5');
lgd = legend('n=0.5',strcat('slope=',num2str(slope(2,1))));
lgd.Location = 'southwest';
%saveas(gcf,'12_L100_rep10000_2tau.png')
hold off

%filename = "12_L100_rep10000_2tau.mat";
%save(filename)