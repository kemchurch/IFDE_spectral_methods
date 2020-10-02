u1 = 10;
N = 32;
load('h2_betastar_proven.mat')
plot(betacont,epscont,'red')
hold on
load('h0_2_h1_proven.mat')
eh_02_h1 = zeros(length(B),1);
for i=1:length(eh_02_h1)
    [a,e] = desingularize_transform(A(:,i),N,u1);
    eh_02_h1(i) = e;
end
plot(B,eh_02_h1,'k-')
load('h1_2_h2_proven.mat')
eh1_2_h2 = zeros(length(B),1);
for i=1:length(eh1_2_h2)
    [a,e] = desingularize_transform(A(:,i),N,u1);
    eh1_2_h2(i) = e;
end
plot(B,eh1_2_h2,'k-')
load('pastbetastar_segment1.mat')
plot(betacont_more,epscont_more,'b--')
load('pastbetastar_segment2.mat')
plot(betacont_more,epscont_more,'b--')