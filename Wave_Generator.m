% Matlab code to generate the random snapshot for Laplace-like pdf:
M = 50;                                                                % Number of snapshots
%N = 20;                      % for Laplace-like                 % Number of large-scale Gaussians
N = 150;                          % for Gaussian-like
Peak = 100;                                                           % Peak of the 3D signal
Sigma = 10;                                                           % std of large scale Gaussian
Ns = 150;                                                            % Number of the small-scale Gaussians
Ss = 3;                     % for Laplace-like                  % std of the small scale Gaussians
%Ss = 10;                        % for Gaussian-like

V = 100 * (0 : 0.01 : 1);
[X,Y] = meshgrid(V);

z = [];
P = rand(N,3);
P(:,1:2) = Peak * P(:,1:2);
for Wnum = 1 : M
    Z = zeros(size(X));
    rand('state', sum(100*clock));
    
    % Generation of the large scale Gaussians
    %P = Peak * rand(N,3);
    for k = 1: N
        Z = Z + 2 * Peak * (P(k,3) - 0.5) / sqrt(2*pi) / Sigma * exp(-((X-P(k,1)).^2 + (Y-P(k,2)).^2)/2/Sigma^2);
    end,
    
    % Genmeration of small scale Gaussians
    Ps = rand(Ns,3);
    Ps(:,1:2) = Peak * Ps(:,1:2);   
    for k = 1 : Ns
        Z = Z + 5 * (Ps(k,3) - 0.5) / sqrt(2*pi) / Ss * exp(-((X-Ps(k,1)).^2 + (Y-Ps(k,2)).^2)/2/Ss^2);
    end,
    
    % Snapshot normalization
    Z = Peak / 2 * Z/(max(max(Z)) - min(min(Z))) + 50;
    
    Data = Z;
    %[p,xx] = ksdensity(Z(:));
    [p,xx] = hist(Z(:),500);
    % saving snapshot
    name = sprintf('NewWave_%d',Wnum);
    save(name,'Data', 'P','Ps','Sigma','Ss','N','Ns','Peak', 'xx', 'p'); 
    z = [z;Z(:)];
end,

%[N,xx] = hist(z,100);
%p = N/sum(N);
%[p,xx] = ksdensity(z);
[p,xx] = hist(Z(:),500);
save('pdf','p','xx');
%figure(1); plot(xx,p); grid on;
%t = 0:100;
%pg=1/sqrt(2*pi)/13*exp(-(t-50).^2/13^2);
%figure(1); hold on; plot(t,pg,'r');

figure(2); mesh(X,Y,Z); hold on;
Vec = 5:10:95;
contour3(X,Y,Z,Vec);


