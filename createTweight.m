function [WT] = createTweight(wBT,AT,MT,n)
% CREATETWEIGHT returns the weight function necessary to impose the desired 
%               bandwidth wBT, the maximum steady-state error AT, the maximum
%               peak MT and the roll-of n 

wT1 = tf([1 wBT],[AT^(1/n) wBT*MT^(1/n)]);
WT = tf(1);

for i=1:n
    WT = WT*wT1;
end
    
end