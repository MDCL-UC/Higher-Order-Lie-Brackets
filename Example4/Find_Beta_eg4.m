clc
clear

params

% Betas
beta2 = beta2_fun(k,p_val, w_val);

beta3 = beta3_fun(k,p_val, w_val);

beta_val = beta_fun(beta2,beta3);

function out =beta_fun(beta2,beta3)

%% Betas
beta = dictionary();
beta(1211) = beta2(1211) - beta3(1211);
beta(1212) = beta2(1212) - beta3(1212);
beta(1213) = beta2(1213) - beta3(1213);
beta(1214) =  beta2(1214) - beta3(1214);

beta(1221) = beta2(1221) - beta3(1221);
beta(1222) = beta2(1222) - beta3(1222);
beta(1223) = beta2(1223) - beta3(1223);
beta(1224) = beta2(1224) - beta3(1224);

beta(1231) = beta2(1231) - beta3(1231);
beta(1232) = beta2(1232) - beta3(1232);
beta(1233) = beta2(1233) - beta3(1233);
beta(1234) = beta2(1234) - beta3(1234);

beta(1241) = beta2(1241) - beta3(1241);
beta(1242) = beta2(1242) - beta3(1242);
beta(1243) = beta2(1243) - beta3(1243);
beta(1244) = beta2(1244) - beta3(1244);

% % % % %
% % % % beta(1311) = beta2(1311) - beta3(1311);
% % % % beta(1312) = beta2(1312) - beta3(1312);
% % % % beta(1313) = beta2(1313) - beta3(1313);
% % % % beta(1314) = beta2(1314) - beta3(1314);
% % % % 
% % % % beta(1321) = beta2(1321) - beta3(1321);
% % % % beta(1322) = beta2(1322) - beta3(1322);
% % % % beta(1323) = beta2(1323) - beta3(1323);
% % % % beta(1324) = beta2(1324) - beta3(1324);
% % % % 
% % % % beta(1331) = beta2(1331) - beta3(1331);
% % % % beta(1332) = beta2(1332) - beta3(1332);
% % % % beta(1333) = beta2(1333) - beta3(1333);
% % % % beta(1334) = beta2(1334) - beta3(1334);
% % % % 
% % % % beta(1341) = beta2(1341) - beta3(1341);
% % % % beta(1342) = beta2(1342) - beta3(1342);
% % % % beta(1343) = beta2(1343) - beta3(1343);
% % % % beta(1344) = beta2(1344) - beta3(1344);
%

beta(1411) = beta2(1411) - beta3(1411);
beta(1412) = beta2(1412) - beta3(1412);
beta(1413) = beta2(1413) - beta3(1413);
beta(1414) = beta2(1414) - beta3(1414);

beta(1421) = beta2(1421) - beta3(1421);
beta(1422) = beta2(1422) - beta3(1422);
beta(1423) = beta2(1423) - beta3(1423);
beta(1424) = beta2(1424) - beta3(1424);

beta(1431) = beta2(1431) - beta3(1431);
beta(1432) = beta2(1432) - beta3(1432);
beta(1433) = beta2(1433) - beta3(1433);
beta(1434) = beta2(1434) - beta3(1434);

beta(1441) = beta2(1441) - beta3(1441);
beta(1442) = beta2(1442) - beta3(1442);
beta(1443) = beta2(1443) - beta3(1443);
beta(1444) = beta2(1444) - beta3(1444);

%
beta(2311) = beta2(2311) - beta3(2311);
beta(2312) = beta2(2312) - beta3(2312);
beta(2313) = beta2(2313) - beta3(2313);
beta(2314) = beta2(2314) - beta3(2314);

beta(2321) = beta2(2321) - beta3(2321);
beta(2322) = beta2(2322) - beta3(2322);
beta(2323) = beta2(2323) - beta3(2323);
beta(2324) = beta2(2324) - beta3(2324);

beta(2331) = beta2(2331) - beta3(2331);
beta(2332) = beta2(2332) - beta3(2332);
beta(2333) = beta2(2333) - beta3(2333);
beta(2334) = beta2(2334) - beta3(2334);

beta(2341) = beta2(2341) - beta3(2341);
beta(2342) = beta2(2342) - beta3(2342);
beta(2343) = beta2(2343) - beta3(2343);
beta(2344) = beta2(2344) - beta3(2344);
%

beta(3411) = beta2(3411) - beta3(3411);
beta(3412) = beta2(3412) - beta3(3412);
beta(3413) = beta2(3413) - beta3(3413);
beta(3414) = beta2(3414) - beta3(3414);

beta(3421) = beta2(3421) - beta3(3421);
beta(3422) = beta2(3422) - beta3(3422);
beta(3423) = beta2(3423) - beta3(3423);
beta(3424) = beta2(3424) - beta3(3424);

beta(3431) = beta2(3431) - beta3(3431);
beta(3432) = beta2(3432) - beta3(3432);
beta(3433) = beta2(3433) - beta3(3433);
beta(3434) = beta2(3434) - beta3(3434);

beta(3441) = beta2(3441) - beta3(3441);
beta(3442) = beta2(3442) - beta3(3442);
beta(3443) = beta2(3443) - beta3(3443);
beta(3444) = beta2(3444) - beta3(3444);

%% Output
out = beta;
beta_val = out;
save("beta_val","beta_val")

end


function out = beta2_fun(k,p_val, w)

%% Parameters
k1 = k(1);
k2 = k(2);
k3 = k(3);
k4 = k(4);

%% Control inputs
syms p q r tau 

% u
u1_r = cos(k1*r);
u2_r = sin(k2*r);
u3_r = cos(k3*r);
u4_r = sin(k4*r);

u1_q = cos(k1*q);
u2_q = sin(k2*q);
u3_q = cos(k3*q);
u4_q = sin(k4*q);

u1_p = cos(k1*p);
u2_p = sin(k2*p);
u3_p = cos(k3*p);
u4_p = sin(k4*p);

u1_tau = cos(k1*tau);
u2_tau = sin(k2*tau);
u3_tau = cos(k3*tau);
u4_tau = sin(k4*tau);


% U
U1_q = int(u1_r,0,q);
U2_q = int(u2_r,0,q);
U3_q = int(u3_r,0,q);
U4_q = int(u4_r,0,q);



%% Alphas
% Alpha1
alpha1_12 =  u2_q*U1_q-u1_q*U2_q;
alpha1_13 =  u3_q*U1_q-u1_q*U3_q;
alpha1_14 =  u4_q*U1_q-u1_q*U4_q;
alpha1_23 =  u3_q*U2_q-u2_q*U3_q;
alpha1_24 =  u4_q*U2_q-u2_q*U4_q;
alpha1_34 =  u4_q*U3_q-u3_q*U4_q;

% Alpha2
alpha2_12 = int(alpha1_12,q,0,p);
alpha2_13 = int(alpha1_13,q,0,p);
alpha2_14 = int(alpha1_14,q,0,p);
alpha2_23 = int(alpha1_23,q,0,p);
alpha2_24 = int(alpha1_24,q,0,p);
alpha2_34 = int(alpha1_34,q,0,p);

% Alpha 6
alpha6_121 = alpha2_12*u1_p;
alpha6_122 = alpha2_12*u2_p;
alpha6_123 = alpha2_12*u3_p;
alpha6_124 = alpha2_12*u4_p;

alpha6_131 = alpha2_13*u1_p;
alpha6_132 = alpha2_13*u2_p;
alpha6_133 = alpha2_13*u3_p;
alpha6_134 = alpha2_13*u4_p;

alpha6_141 = alpha2_14*u1_p;
alpha6_142 = alpha2_14*u2_p;
alpha6_143 = alpha2_14*u3_p;
alpha6_144 = alpha2_14*u4_p;

alpha6_231 = alpha2_23*u1_p;
alpha6_232 = alpha2_23*u2_p;
alpha6_233 = alpha2_23*u3_p;
alpha6_234 = alpha2_23*u4_p;

alpha6_341 = alpha2_34*u1_p;
alpha6_342 = alpha2_34*u2_p;
alpha6_343 = alpha2_34*u3_p;
alpha6_344 = alpha2_34*u4_p;

% Alpha 7 
alpha7_121 = int(alpha6_121,p,0,tau);
alpha7_122 = int(alpha6_122,p,0,tau);
alpha7_123 = int(alpha6_123,p,0,tau);
alpha7_124 = int(alpha6_124,p,0,tau);

alpha7_131 = int(alpha6_131,p,0,tau);
alpha7_132 = int(alpha6_132,p,0,tau);
alpha7_133 = int(alpha6_133,p,0,tau);
alpha7_134 = int(alpha6_134,p,0,tau);

alpha7_141 = int(alpha6_141,p,0,tau);
alpha7_142 = int(alpha6_142,p,0,tau);
alpha7_143 = int(alpha6_143,p,0,tau);
alpha7_144 = int(alpha6_144,p,0,tau);

alpha7_231 = int(alpha6_231,p,0,tau);
alpha7_232 = int(alpha6_232,p,0,tau);
alpha7_233 = int(alpha6_233,p,0,tau);
alpha7_234 = int(alpha6_234,p,0,tau);

alpha7_341 = int(alpha6_341,p,0,tau);
alpha7_342 = int(alpha6_342,p,0,tau);
alpha7_343 = int(alpha6_343,p,0,tau);
alpha7_344 = int(alpha6_344,p,0,tau);

%% Alpha 8
% Alpha 8
% 12
alpha8_1211 = alpha7_121*u1_tau;
alpha8_1212 = alpha7_121*u2_tau;
alpha8_1213 = alpha7_121*u3_tau;
alpha8_1214 = alpha7_121*u4_tau;

alpha8_1221 = alpha7_122*u1_tau;
alpha8_1222 = alpha7_122*u2_tau;
alpha8_1223 = alpha7_122*u3_tau;
alpha8_1224 = alpha7_122*u4_tau;

alpha8_1231 = alpha7_123*u1_tau;
alpha8_1232 = alpha7_123*u2_tau;
alpha8_1233 = alpha7_123*u3_tau;
alpha8_1234 = alpha7_123*u4_tau;

alpha8_1241 = alpha7_124*u1_tau;
alpha8_1242 = alpha7_124*u2_tau;
alpha8_1243 = alpha7_124*u3_tau;
alpha8_1244 = alpha7_124*u4_tau;

% 13
alpha8_1311 = alpha7_131*u1_tau;
alpha8_1312 = alpha7_131*u2_tau;
alpha8_1313 = alpha7_131*u3_tau;
alpha8_1314 = alpha7_131*u4_tau;


alpha8_1321 = alpha7_132*u1_tau;
alpha8_1322 = alpha7_132*u2_tau;
alpha8_1323 = alpha7_132*u3_tau;
alpha8_1324 = alpha7_132*u4_tau;


alpha8_1331 = alpha7_133*u1_tau;
alpha8_1332 = alpha7_133*u2_tau;
alpha8_1333 = alpha7_133*u3_tau;
alpha8_1334 = alpha7_133*u4_tau;

alpha8_1341 = alpha7_134*u1_tau;
alpha8_1342 = alpha7_134*u2_tau;
alpha8_1343 = alpha7_134*u3_tau;
alpha8_1344 = alpha7_134*u4_tau;

% 14
alpha8_1411 = alpha7_141*u1_tau;
alpha8_1412 = alpha7_141*u2_tau;
alpha8_1413 = alpha7_141*u3_tau;
alpha8_1414 = alpha7_141*u4_tau;


alpha8_1421 = alpha7_142*u1_tau;
alpha8_1422 = alpha7_142*u2_tau;
alpha8_1423 = alpha7_142*u3_tau;
alpha8_1424 = alpha7_142*u4_tau;


alpha8_1431 = alpha7_143*u1_tau;
alpha8_1432 = alpha7_143*u2_tau;
alpha8_1433 = alpha7_143*u3_tau;
alpha8_1434 = alpha7_143*u4_tau;

alpha8_1441 = alpha7_144*u1_tau;
alpha8_1442 = alpha7_144*u2_tau;
alpha8_1443 = alpha7_144*u3_tau;
alpha8_1444 = alpha7_144*u4_tau;


% 23
alpha8_2311 = alpha7_231*u1_tau;
alpha8_2312 = alpha7_231*u2_tau;
alpha8_2313 = alpha7_231*u3_tau;
alpha8_2314 = alpha7_231*u4_tau;


alpha8_2321 = alpha7_232*u1_tau;
alpha8_2322 = alpha7_232*u2_tau;
alpha8_2323 = alpha7_232*u3_tau;
alpha8_2324 = alpha7_232*u4_tau;


alpha8_2331 = alpha7_233*u1_tau;
alpha8_2332 = alpha7_233*u2_tau;
alpha8_2333 = alpha7_233*u3_tau;
alpha8_2334 = alpha7_233*u4_tau;

alpha8_2341 = alpha7_234*u1_tau;
alpha8_2342 = alpha7_234*u2_tau;
alpha8_2343 = alpha7_234*u3_tau;
alpha8_2344 = alpha7_234*u4_tau;

% 34
alpha8_3411 = alpha7_341*u1_tau;
alpha8_3412 = alpha7_341*u2_tau;
alpha8_3413 = alpha7_341*u3_tau;
alpha8_3414 = alpha7_341*u4_tau;


alpha8_3421 = alpha7_342*u1_tau;
alpha8_3422 = alpha7_342*u2_tau;
alpha8_3423 = alpha7_342*u3_tau;
alpha8_3424 = alpha7_342*u4_tau;


alpha8_3431 = alpha7_343*u1_tau;
alpha8_3432 = alpha7_343*u2_tau;
alpha8_3433 = alpha7_343*u3_tau;
alpha8_3434 = alpha7_343*u4_tau;

alpha8_3441 = alpha7_344*u1_tau;
alpha8_3442 = alpha7_344*u2_tau;
alpha8_3443 = alpha7_344*u3_tau;
alpha8_3444 = alpha7_344*u4_tau;


%% Betas
%12
beta2 = dictionary();
beta2(1211) = calculate_beta(w,k,p_val,[1,2,1,1],alpha8_1211);
beta2(1212) = calculate_beta(w,k,p_val,[1,2,1,2],alpha8_1212);
beta2(1213) = calculate_beta(w,k,p_val,[1,2,1,3],alpha8_1213);
beta2(1214) = calculate_beta(w,k,p_val,[1,2,1,4],alpha8_1214);

beta2(1221) = calculate_beta(w,k,p_val,[1,2,2,1],alpha8_1221);
beta2(1222) = calculate_beta(w,k,p_val,[1,2,2,2],alpha8_1222);
beta2(1223) = calculate_beta(w,k,p_val,[1,2,2,3],alpha8_1223);
beta2(1224) = calculate_beta(w,k,p_val,[1,2,2,4],alpha8_1224);

beta2(1231) = calculate_beta(w,k,p_val,[1,2,3,1],alpha8_1231);
beta2(1232) = calculate_beta(w,k,p_val,[1,2,3,2],alpha8_1232);
beta2(1233) = calculate_beta(w,k,p_val,[1,2,3,3],alpha8_1233);
beta2(1234) = calculate_beta(w,k,p_val,[1,2,3,4],alpha8_1234);

beta2(1241) = calculate_beta(w,k,p_val,[1,2,4,1],alpha8_1241);
beta2(1242) = calculate_beta(w,k,p_val,[1,2,4,2],alpha8_1242);
beta2(1243) = calculate_beta(w,k,p_val,[1,2,4,3],alpha8_1243);
beta2(1244) = calculate_beta(w,k,p_val,[1,2,4,4],alpha8_1244);


%14
beta2(1411) = calculate_beta(w,k,p_val,[1,4,1,1],alpha8_1411);
beta2(1412) = calculate_beta(w,k,p_val,[1,4,1,2],alpha8_1412);
beta2(1413) = calculate_beta(w,k,p_val,[1,4,1,3],alpha8_1413);
beta2(1414) = calculate_beta(w,k,p_val,[1,4,1,4],alpha8_1414);

beta2(1421) = calculate_beta(w,k,p_val,[1,4,2,1],alpha8_1421);
beta2(1422) = calculate_beta(w,k,p_val,[1,4,2,2],alpha8_1422);
beta2(1423) = calculate_beta(w,k,p_val,[1,4,2,3],alpha8_1423);
beta2(1424) = calculate_beta(w,k,p_val,[1,4,2,4],alpha8_1424);

beta2(1431) = calculate_beta(w,k,p_val,[1,4,3,1],alpha8_1431);
beta2(1432) = calculate_beta(w,k,p_val,[1,4,3,2],alpha8_1432);
beta2(1433) = calculate_beta(w,k,p_val,[1,4,3,3],alpha8_1433);
beta2(1434) = calculate_beta(w,k,p_val,[1,4,3,4],alpha8_1434);

beta2(1441) = calculate_beta(w,k,p_val,[1,4,4,1],alpha8_1441);
beta2(1442) = calculate_beta(w,k,p_val,[1,4,4,2],alpha8_1442);
beta2(1443) = calculate_beta(w,k,p_val,[1,4,4,3],alpha8_1443);
beta2(1444) = calculate_beta(w,k,p_val,[1,4,4,4],alpha8_1444);

%23
beta2(2311) = calculate_beta(w,k,p_val,[2,3,1,1],alpha8_2311);
beta2(2312) = calculate_beta(w,k,p_val,[2,3,1,2],alpha8_2312);
beta2(2313) = calculate_beta(w,k,p_val,[2,3,1,3],alpha8_2313);
beta2(2314) = calculate_beta(w,k,p_val,[2,3,1,4],alpha8_2314);

beta2(2321) = calculate_beta(w,k,p_val,[2,3,2,1],alpha8_2321);
beta2(2322) = calculate_beta(w,k,p_val,[2,3,2,2],alpha8_2322);
beta2(2323) = calculate_beta(w,k,p_val,[2,3,2,3],alpha8_2323);
beta2(2324) = calculate_beta(w,k,p_val,[2,3,2,4],alpha8_2324);

beta2(2331) = calculate_beta(w,k,p_val,[2,3,3,1],alpha8_2331);
beta2(2332) = calculate_beta(w,k,p_val,[2,3,3,2],alpha8_2332);
beta2(2333) = calculate_beta(w,k,p_val,[2,3,3,3],alpha8_2333);
beta2(2334) = calculate_beta(w,k,p_val,[2,3,3,4],alpha8_2334);

beta2(2341) = calculate_beta(w,k,p_val,[2,3,4,1],alpha8_2341);
beta2(2342) = calculate_beta(w,k,p_val,[2,3,4,2],alpha8_2342);
beta2(2343) = calculate_beta(w,k,p_val,[2,3,4,3],alpha8_2343);
beta2(2344) = calculate_beta(w,k,p_val,[2,3,4,4],alpha8_2344);

%34
beta2(3411) = calculate_beta(w,k,p_val,[3,4,1,1],alpha8_3411);
beta2(3412) = calculate_beta(w,k,p_val,[3,4,1,2],alpha8_3412);
beta2(3413) = calculate_beta(w,k,p_val,[3,4,1,3],alpha8_3413);
beta2(3414) = calculate_beta(w,k,p_val,[3,4,1,4],alpha8_3414);

beta2(3421) = calculate_beta(w,k,p_val,[3,4,2,1],alpha8_3421);
beta2(3422) = calculate_beta(w,k,p_val,[3,4,2,2],alpha8_3422);
beta2(3423) = calculate_beta(w,k,p_val,[3,4,2,3],alpha8_3423);
beta2(3424) = calculate_beta(w,k,p_val,[3,4,2,4],alpha8_3424);

beta2(3431) = calculate_beta(w,k,p_val,[3,4,3,1],alpha8_3431);
beta2(3432) = calculate_beta(w,k,p_val,[3,4,3,2],alpha8_3432);
beta2(3433) = calculate_beta(w,k,p_val,[3,4,3,3],alpha8_3433);
beta2(3434) = calculate_beta(w,k,p_val,[3,4,3,4],alpha8_3434);

beta2(3441) = calculate_beta(w,k,p_val,[3,4,4,1],alpha8_3441);
beta2(3442) = calculate_beta(w,k,p_val,[3,4,4,2],alpha8_3442);
beta2(3443) = calculate_beta(w,k,p_val,[3,4,4,3],alpha8_3443);
beta2(3444) = calculate_beta(w,k,p_val,[3,4,4,4],alpha8_3444);

%% Output
out = beta2;
save("beta2","beta2");

end


function out = beta3_fun(k,p_val, w)

k1 = k(1);
k2 = k(2);
k3 = k(3);
k4 = k(4);

%% Control inputs
syms p q r tau 

u1_r = cos(k1*r);
u2_r = sin(k2*r);
u3_r = cos(k3*r);
u4_r = sin(k4*r);

u1_q = cos(k1*q);
u2_q = sin(k2*q);
u3_q = cos(k3*q);
u4_q = sin(k4*q);

u1_p = cos(k1*p);
u2_p = sin(k2*p);
u3_p = cos(k3*p);
u4_p = sin(k4*p);

u1_tau = cos(k1*tau);
u2_tau = sin(k2*tau);
u3_tau = cos(k3*tau);
u4_tau = sin(k4*tau);

% U
U1_p = int(u1_q,0,p);
U2_p = int(u2_q,0,p);
U3_p = int(u3_q,0,p);
U4_p = int(u4_q,0,p);
%% Alphas
% Alpha 1
alpha1_12 =  u2_p*U1_p-u1_p*U2_p;
alpha1_13 =  u3_p*U1_p-u1_p*U3_p;
alpha1_14 =  u4_p*U1_p-u1_p*U4_p;
alpha1_23 =  u3_p*U2_p-u2_p*U3_p;
alpha1_24 =  u4_p*U2_p-u2_p*U4_p;
alpha1_34 =  u4_p*U3_p-u3_p*U4_p;


% Alpha 9
alpha9_121 = alpha1_12*u1_tau;
alpha9_122 = alpha1_12*u2_tau;
alpha9_123 = alpha1_12*u3_tau;
alpha9_124 = alpha1_12*u4_tau;

alpha9_131 = alpha1_13*u1_tau;
alpha9_132 = alpha1_13*u2_tau;
alpha9_133 = alpha1_13*u3_tau;
alpha9_134 = alpha1_13*u4_tau;

alpha9_141 = alpha1_14*u1_tau;
alpha9_142 = alpha1_14*u2_tau;
alpha9_143 = alpha1_14*u3_tau;
alpha9_144 = alpha1_14*u4_tau;

alpha9_231 = alpha1_23*u1_tau;
alpha9_232 = alpha1_23*u2_tau;
alpha9_233 = alpha1_23*u3_tau;
alpha9_234 = alpha1_23*u4_tau;

alpha9_341 = alpha1_34*u1_tau;
alpha9_342 = alpha1_34*u2_tau;
alpha9_343 = alpha1_34*u3_tau;
alpha9_344 = alpha1_34*u4_tau;

% Alpha 10
alpha10_1211 = int(alpha9_121*U1_p,p,0,tau);
alpha10_1212 = int(alpha9_121*U2_p,p,0,tau);
alpha10_1213 = int(alpha9_121*U3_p,p,0,tau);
alpha10_1214 = int(alpha9_121*U4_p,p,0,tau);

alpha10_1221 = int(alpha9_122*U1_p,p,0,tau);
alpha10_1222 = int(alpha9_122*U2_p,p,0,tau);
alpha10_1223 = int(alpha9_122*U3_p,p,0,tau);
alpha10_1224 = int(alpha9_122*U4_p,p,0,tau);

alpha10_1231 = int(alpha9_123*U1_p,p,0,tau);
alpha10_1232 = int(alpha9_123*U2_p,p,0,tau);
alpha10_1233 = int(alpha9_123*U3_p,p,0,tau);
alpha10_1234 = int(alpha9_123*U4_p,p,0,tau);

alpha10_1241 = int(alpha9_124*U1_p,p,0,tau);
alpha10_1242 = int(alpha9_124*U2_p,p,0,tau);
alpha10_1243 = int(alpha9_124*U3_p,p,0,tau);
alpha10_1244 = int(alpha9_124*U4_p,p,0,tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha10_1311 = int(alpha9_131*U1_p,p,0,tau);
alpha10_1312 = int(alpha9_131*U2_p,p,0,tau);
alpha10_1313 = int(alpha9_131*U3_p,p,0,tau);
alpha10_1314 = int(alpha9_131*U4_p,p,0,tau);

alpha10_1321 = int(alpha9_132*U1_p,p,0,tau);
alpha10_1322 = int(alpha9_132*U2_p,p,0,tau);
alpha10_1323 = int(alpha9_132*U3_p,p,0,tau);
alpha10_1324 = int(alpha9_132*U4_p,p,0,tau);

alpha10_1331 = int(alpha9_133*U1_p,p,0,tau);
alpha10_1332 = int(alpha9_133*U2_p,p,0,tau);
alpha10_1333 = int(alpha9_133*U3_p,p,0,tau);
alpha10_1334 = int(alpha9_133*U4_p,p,0,tau);

alpha10_1341 = int(alpha9_134*U1_p,p,0,tau);
alpha10_1342 = int(alpha9_134*U2_p,p,0,tau);
alpha10_1343 = int(alpha9_134*U3_p,p,0,tau);
alpha10_1344 = int(alpha9_134*U4_p,p,0,tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha10_1411 = int(alpha9_141*U1_p,p,0,tau);
alpha10_1412 = int(alpha9_141*U2_p,p,0,tau);
alpha10_1413 = int(alpha9_141*U3_p,p,0,tau);
alpha10_1414 = int(alpha9_141*U4_p,p,0,tau);

alpha10_1421 = int(alpha9_142*U1_p,p,0,tau);
alpha10_1422 = int(alpha9_142*U2_p,p,0,tau);
alpha10_1423 = int(alpha9_142*U3_p,p,0,tau);
alpha10_1424 = int(alpha9_142*U4_p,p,0,tau);

alpha10_1431 = int(alpha9_143*U1_p,p,0,tau);
alpha10_1432 = int(alpha9_143*U2_p,p,0,tau);
alpha10_1433 = int(alpha9_143*U3_p,p,0,tau);
alpha10_1434 = int(alpha9_143*U4_p,p,0,tau);

alpha10_1441 = int(alpha9_144*U1_p,p,0,tau);
alpha10_1442 = int(alpha9_144*U2_p,p,0,tau);
alpha10_1443 = int(alpha9_144*U3_p,p,0,tau);
alpha10_1444 = int(alpha9_144*U4_p,p,0,tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha10_2311 = int(alpha9_231*U1_p,p,0,tau);
alpha10_2312 = int(alpha9_231*U2_p,p,0,tau);
alpha10_2313 = int(alpha9_231*U3_p,p,0,tau);
alpha10_2314 = int(alpha9_231*U4_p,p,0,tau);

alpha10_2321 = int(alpha9_232*U1_p,p,0,tau);
alpha10_2322 = int(alpha9_232*U2_p,p,0,tau);
alpha10_2323 = int(alpha9_232*U3_p,p,0,tau);
alpha10_2324 = int(alpha9_232*U4_p,p,0,tau);

alpha10_2331 = int(alpha9_233*U1_p,p,0,tau);
alpha10_2332 = int(alpha9_233*U2_p,p,0,tau);
alpha10_2333 = int(alpha9_233*U3_p,p,0,tau);
alpha10_2334 = int(alpha9_233*U4_p,p,0,tau);

alpha10_2341 = int(alpha9_234*U1_p,p,0,tau);
alpha10_2342 = int(alpha9_234*U2_p,p,0,tau);
alpha10_2343 = int(alpha9_234*U3_p,p,0,tau);
alpha10_2344 = int(alpha9_234*U4_p,p,0,tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha10_3411 = int(alpha9_341*U1_p,p,0,tau);
alpha10_3412 = int(alpha9_341*U2_p,p,0,tau);
alpha10_3413 = int(alpha9_341*U3_p,p,0,tau);
alpha10_3414 = int(alpha9_341*U4_p,p,0,tau);

alpha10_3421 = int(alpha9_342*U1_p,p,0,tau);
alpha10_3422 = int(alpha9_342*U2_p,p,0,tau);
alpha10_3423 = int(alpha9_342*U3_p,p,0,tau);
alpha10_3424 = int(alpha9_342*U4_p,p,0,tau);

alpha10_3431 = int(alpha9_343*U1_p,p,0,tau);
alpha10_3432 = int(alpha9_343*U2_p,p,0,tau);
alpha10_3433 = int(alpha9_343*U3_p,p,0,tau);
alpha10_3434 = int(alpha9_343*U4_p,p,0,tau);

alpha10_3441 = int(alpha9_344*U1_p,p,0,tau);
alpha10_3442 = int(alpha9_344*U2_p,p,0,tau);
alpha10_3443 = int(alpha9_344*U3_p,p,0,tau);
alpha10_3444 = int(alpha9_344*U4_p,p,0,tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Betas
%12
beta3 = dictionary();
beta3(1211) = calculate_beta(w,k,p_val,[1,2,1,1],alpha10_1211);
beta3(1212) = calculate_beta(w,k,p_val,[1,2,1,2],alpha10_1212);
beta3(1213) = calculate_beta(w,k,p_val,[1,2,1,3],alpha10_1213);
beta3(1214) = calculate_beta(w,k,p_val,[1,2,1,4],alpha10_1214);

beta3(1221) = calculate_beta(w,k,p_val,[1,2,2,1],alpha10_1221);
beta3(1222) = calculate_beta(w,k,p_val,[1,2,2,2],alpha10_1222);
beta3(1223) = calculate_beta(w,k,p_val,[1,2,2,3],alpha10_1223);
beta3(1224) = calculate_beta(w,k,p_val,[1,2,2,4],alpha10_1224);

beta3(1231) = calculate_beta(w,k,p_val,[1,2,3,1],alpha10_1231);
beta3(1232) = calculate_beta(w,k,p_val,[1,2,3,2],alpha10_1232);
beta3(1233) = calculate_beta(w,k,p_val,[1,2,3,3],alpha10_1233);
beta3(1234) = calculate_beta(w,k,p_val,[1,2,3,4],alpha10_1234);

beta3(1241) = calculate_beta(w,k,p_val,[1,2,4,1],alpha10_1241);
beta3(1242) = calculate_beta(w,k,p_val,[1,2,4,2],alpha10_1242);
beta3(1243) = calculate_beta(w,k,p_val,[1,2,4,3],alpha10_1243);
beta3(1244) = calculate_beta(w,k,p_val,[1,2,4,4],alpha10_1244);


% % % % %13
% % % % beta3(1311) = calculate_beta(w,k,p_val,[1,3,1,1],alpha10_1311);
% % % % beta3(1312) = calculate_beta(w,k,p_val,[1,3,1,2],alpha10_1312);
% % % % beta3(1313) = calculate_beta(w,k,p_val,[1,3,1,3],alpha10_1313);
% % % % beta3(1314) = calculate_beta(w,k,p_val,[1,3,1,4],alpha10_1314);
% % % % 
% % % % beta3(1321) = calculate_beta(w,k,p_val,[1,3,2,1],alpha10_1321);
% % % % beta3(1322) = calculate_beta(w,k,p_val,[1,3,2,2],alpha10_1322);
% % % % beta3(1323) = calculate_beta(w,k,p_val,[1,3,2,3],alpha10_1323);
% % % % beta3(1324) = calculate_beta(w,k,p_val,[1,3,2,4],alpha10_1324);
% % % % 
% % % % beta3(1331) = calculate_beta(w,k,p_val,[1,3,3,1],alpha10_1331);
% % % % beta3(1332) = calculate_beta(w,k,p_val,[1,3,3,2],alpha10_1332);
% % % % beta3(1333) = calculate_beta(w,k,p_val,[1,3,3,3],alpha10_1333);
% % % % beta3(1334) = calculate_beta(w,k,p_val,[1,3,3,4],alpha10_1334);
% % % % 
% % % % beta3(1341) = calculate_beta(w,k,p_val,[1,3,4,1],alpha10_1341);
% % % % beta3(1342) = calculate_beta(w,k,p_val,[1,3,4,2],alpha10_1342);
% % % % beta3(1343) = calculate_beta(w,k,p_val,[1,3,4,3],alpha10_1343);
% % % % beta3(1344) = calculate_beta(w,k,p_val,[1,3,4,4],alpha10_1344);

%14
beta3(1411) = calculate_beta(w,k,p_val,[1,4,1,1],alpha10_1411);
beta3(1412) = calculate_beta(w,k,p_val,[1,4,1,2],alpha10_1412);
beta3(1413) = calculate_beta(w,k,p_val,[1,4,1,3],alpha10_1413);
beta3(1414) = calculate_beta(w,k,p_val,[1,4,1,4],alpha10_1414);

beta3(1421) = calculate_beta(w,k,p_val,[1,4,2,1],alpha10_1421);
beta3(1422) = calculate_beta(w,k,p_val,[1,4,2,2],alpha10_1422);
beta3(1423) = calculate_beta(w,k,p_val,[1,4,2,3],alpha10_1423);
beta3(1424) = calculate_beta(w,k,p_val,[1,4,2,4],alpha10_1424);

beta3(1431) = calculate_beta(w,k,p_val,[1,4,3,1],alpha10_1431);
beta3(1432) = calculate_beta(w,k,p_val,[1,4,3,2],alpha10_1432);
beta3(1433) = calculate_beta(w,k,p_val,[1,4,3,3],alpha10_1433);
beta3(1434) = calculate_beta(w,k,p_val,[1,4,3,4],alpha10_1434);

beta3(1441) = calculate_beta(w,k,p_val,[1,4,4,1],alpha10_1441);
beta3(1442) = calculate_beta(w,k,p_val,[1,4,4,2],alpha10_1442);
beta3(1443) = calculate_beta(w,k,p_val,[1,4,4,3],alpha10_1443);
beta3(1444) = calculate_beta(w,k,p_val,[1,4,4,4],alpha10_1444);

%23
beta3(2311) = calculate_beta(w,k,p_val,[2,3,1,1],alpha10_2311);
beta3(2312) = calculate_beta(w,k,p_val,[2,3,1,2],alpha10_2312);
beta3(2313) = calculate_beta(w,k,p_val,[2,3,1,3],alpha10_2313);
beta3(2314) = calculate_beta(w,k,p_val,[2,3,1,4],alpha10_2314);

beta3(2321) = calculate_beta(w,k,p_val,[2,3,2,1],alpha10_2321);
beta3(2322) = calculate_beta(w,k,p_val,[2,3,2,2],alpha10_2322);
beta3(2323) = calculate_beta(w,k,p_val,[2,3,2,3],alpha10_2323);
beta3(2324) = calculate_beta(w,k,p_val,[2,3,2,4],alpha10_2324);

beta3(2331) = calculate_beta(w,k,p_val,[2,3,3,1],alpha10_2331);
beta3(2332) = calculate_beta(w,k,p_val,[2,3,3,2],alpha10_2332);
beta3(2333) = calculate_beta(w,k,p_val,[2,3,3,3],alpha10_2333);
beta3(2334) = calculate_beta(w,k,p_val,[2,3,3,4],alpha10_2334);

beta3(2341) = calculate_beta(w,k,p_val,[2,3,4,1],alpha10_2341);
beta3(2342) = calculate_beta(w,k,p_val,[2,3,4,2],alpha10_2342);
beta3(2343) = calculate_beta(w,k,p_val,[2,3,4,3],alpha10_2343);
beta3(2344) = calculate_beta(w,k,p_val,[2,3,4,4],alpha10_2344);

%34
beta3(3411) = calculate_beta(w,k,p_val,[3,4,1,1],alpha10_3411);
beta3(3412) = calculate_beta(w,k,p_val,[3,4,1,2],alpha10_3412);
beta3(3413) = calculate_beta(w,k,p_val,[3,4,1,3],alpha10_3413);
beta3(3414) = calculate_beta(w,k,p_val,[3,4,1,4],alpha10_3414);

beta3(3421) = calculate_beta(w,k,p_val,[3,4,2,1],alpha10_3421);
beta3(3422) = calculate_beta(w,k,p_val,[3,4,2,2],alpha10_3422);
beta3(3423) = calculate_beta(w,k,p_val,[3,4,2,3],alpha10_3423);
beta3(3424) = calculate_beta(w,k,p_val,[3,4,2,4],alpha10_3424);

beta3(3431) = calculate_beta(w,k,p_val,[3,4,3,1],alpha10_3431);
beta3(3432) = calculate_beta(w,k,p_val,[3,4,3,2],alpha10_3432);
beta3(3433) = calculate_beta(w,k,p_val,[3,4,3,3],alpha10_3433);
beta3(3434) = calculate_beta(w,k,p_val,[3,4,3,4],alpha10_3434);

beta3(3441) = calculate_beta(w,k,p_val,[3,4,4,1],alpha10_3441);
beta3(3442) = calculate_beta(w,k,p_val,[3,4,4,2],alpha10_3442);
beta3(3443) = calculate_beta(w,k,p_val,[3,4,4,3],alpha10_3443);
beta3(3444) = calculate_beta(w,k,p_val,[3,4,4,4],alpha10_3444);

%% Output
out = beta3;
save("beta3","beta3")

end



function beta_out = calculate_beta(w,k,p,order,alpha)
syms tau 
m = calculate_m(k,order);
sum_p = calculate_sump(p,order);
T = 2*pi*m;
beta_out = w^sum_p/(12*T)*int(alpha,tau,0,T);
end

function m = calculate_m(k,order)
n = length(k);
switch n
    case 2
       k1 = k(order(1));
       k2 = k(order(2));
       term1 = ceil(1/k1);
       term2 = ceil(1/k2);
       m= lcm(term1,term2);
       
    case 3
       k1 = k(order(1));
       k2 = k(order(2));
       k3 = k(order(3));
       term1 = ceil(1/k1);
       term2 = ceil(1/k2);
       term3 = ceil(1/k3);
       lcm12 = lcm(term1,term2);
       m = lcm(lcm12,term3);

    case 4
       k1 = k(order(1));
       k2 = k(order(2));
       k3 = k(order(3));
       k4 = k(order(4));
       term1 = ceil(1/k1);
       term2 = ceil(1/k2);
       term3 = ceil(1/k3);
       term4 = ceil(1/k4);
       lcm12 = lcm(term1,term2);
       lcm123  = lcm(lcm12,term3);
       m = lcm(lcm123,term4);
end
end

function sump = calculate_sump(p,order)
n = length(order); % Number of elements in p

switch n
    case 2
        i = order(1);
        j = order(2);
        sump = p(i) + p(j) - 1;
    case 3
        i = order(1);
        j = order(2);
        k = order(3);
        sump = p(i) + p(j) + p(k) - 2;

    case 4
        i = order(1);
        j = order(2);
        k = order(3);
        l = order(4);
        sump = p(i) + p(j) + p(k) +p(l) - 3;
    otherwise
        error('Unsupported number of elements in p or incorrect number of indices.');
end
end



