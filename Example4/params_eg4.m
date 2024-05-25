%% Parameters
k1 = 1;           
k2 = 1;           
k3 = 1/4; 
k4 = 1/4;  
k = [k1,k2,k3,k4];


p1 = 0.99 ;
p2 = 0.01 ;
p3 = 0.99 ;
p4 = 0.01;
p_val = [p1,p2,p3,p4]; % Write in order

w_val = 100;

u1 = @(x) cos(x);
u2 = @(x) sin(x);
u3 = @(x) cos(x);
u4 = @(x) sin(x);

% Indices for nu_ij
index_all2 = [12,14,23,34];

% Indices for nu_ijk
index_all3 = [121,122,123,124,...
    141,142,143,144,...
    231,232,233,234,...
    341,342,343,344];
