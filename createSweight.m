function [WS] = createSweight(wB,A,M,n)
% CREATESWEIGHT returns the weight function necessary to impose the desired 
%               bandwidth wB, the maximum steady-state error A, the maximum
%               peak M and the order n (1 for step 2 - for ramp)

wS1 = tf([1/M^(1/n) wB],[1 wB*A^(1/n)]);
WS = tf(1);

for i=1:n
    WS = WS*wS1;
end
    
end

