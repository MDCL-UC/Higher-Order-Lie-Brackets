
%%
function Find_nu


syms p w tau
params


u1_p = u1(k1 * p);      
u2_p = u2(k2 * p);      
u3_p = u3(k3 * p);
u4_p = u4(k4 * p);
u_p = [u1_p, u2_p, u3_p, u4_p];

u1_tau = u1(k1 * tau);     
u2_tau = u2(k2 * tau);    
u3_tau = u3(k3 * tau); 
u4_tau = u4(k4 * tau); 
u_tau = [u1_tau, u2_tau, u3_tau, u4_tau];


% Find the integral of u_p
U1 = int(u1_p, p, 0, tau);
U2 = int(u2_p, p, 0, tau);
U3 = int(u3_p, p, 0, tau);
U4 = int(u4_p, p, 0, tau);
U = [U1, U2, U3, U4];


% For nu_ij
nu2 = dictionary();
for i=1:length(index_all2)
    index = index_all2(i);    
    nu2(index) = find_nu_val(index,k,p_val,u_tau,U,tau,w,w_val);
end
% nu_val2 = double(subs(values(nu2),w,w_val));
nu_val2 = nu2;
save("nu_val2","nu_val2")


% For nu_ijk
nu3 = dictionary();
for i=1:length(index_all3)
    index = index_all3(i);    
    nu3(index) = find_nu_val(index,k,p_val,u_tau,U,tau,w,w_val);
end
% nu_val3 = double(subs(values(nu3),w,w_val));
nu_val3 = nu3;
save("nu_val3","nu_val3")

end



function nu = find_nu_val(index,k,p_val,u_tau,U,tau,w,w_val)
switch 1
    case index>10 & index<100 % 2 digit
        k2_pos = rem(index,10);
        k1_pos = (index-k2_pos)/10;

        k1 = k(k1_pos);
        k2 = k(k2_pos);
        U1 = U(k1_pos);
        U2 = U(k2_pos);
        u1_tau = u_tau(k1_pos);
        u2_tau = u_tau(k2_pos);
        int_val = calculate_int_second_order(k1, k2, U1, U2, u1_tau, u2_tau, tau);
        nu = w_val^calculate_sump(p_val,[k1_pos,k2_pos])*int_val;

    case index>100 && index<1000 % 3 digit
        k3_pos = rem(index,10);
        first2terms = (index-k3_pos)/10;
        k2_pos = rem(first2terms,10);
        k1_pos = (first2terms-k2_pos)/10;

        k1 = k(k1_pos);
        k2 = k(k2_pos);
        k3 = k(k3_pos);
        U1 = U(k1_pos);
        U2 = U(k2_pos);
        U3 = U(k3_pos);
        u1_tau = u_tau(k1_pos);
        u2_tau = u_tau(k2_pos);
        u3_tau = u_tau(k3_pos);
        int_val = calculate_int_third_order(k1, k2, k3, U1, U2, U3, u1_tau, u2_tau, u3_tau, tau);
        nu = w_val^calculate_sump(p_val,[k1_pos,k2_pos,k3_pos])*int_val;
end
end



%%  Necessary functions
function integral = calculate_int_second_order(k1, k2, U1, U2, u1_tau, u2_tau, tau)
    if k1 <= 1 && k2 <= 1
        m = lcm(1 / k1, 1 / k2);
    else
        m = lcm(ceil(1 / k1), ceil(1 / k2));
    end
    T = 2 * pi * m;
    beta = (u2_tau * U1) - (u1_tau * U2);
    integral = eval(1 / (2 * T) * int(beta, tau, 0, T));
end

function integral = calculate_int_third_order(k1, k2,k3, U1, U2, U3, u1_tau, u2_tau,u3_tau,tau)
    if k1<=1 && k2<=1 && k3<=1
        m = lcm(1/k1,lcm(1/k2,1/k3));
    else
        m=lcm(ceil(1/k1),lcm(ceil(1/k2),ceil(1/k3)));
    end
    T = 2 * pi * m;
    beta = (u2_tau * U1) - (u1_tau * U2);
    integral = eval(1 / (3 * T) * int(beta*U3, tau, 0, T));
end


function sump = calculate_sump(p,m)
n = length(m); % Number of elements in p

switch n
    case 2
        i = m(1);
        j = m(2);
        sump = p(i) + p(j) - 1;
    case 3
        i = m(1);
        j = m(2);
        k = m(3);
        sump = p(i) + p(j) + p(k) - 2;

     case 4
        i = m(1);
        j = m(2);
        k = m(3);
        l = m(4);
        sump = p(i) + p(j) + p(k) +p(l) - 3;
    otherwise
        error('Unsupported number of elements in p or incorrect number of indices.');
end
end