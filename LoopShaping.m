
%% Plant model

numG = [0,0,292.423];
denG = [1,14.117,1.259];
G = tf(numG, denG);

 clf, sigma(G,'g',{.1,100});
title('Singular value plot for theeta/U')
s = zpk('s'); 
Gd =10/s;     
sigma(Gd,{.1 100})
grid
title('Target loop shape Gd(s).')
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times')
%%
% create textarrow pointing to crossover frequency Wc
hold on;
plot([10,35],[0,21],'k.-'); 
plot(10,0,'kd');
plot([.1,100],[0 0],'k');
text(3,23,'Crossover Frequency \omega_c = 10');
hold off;
%Using LOOPSYN to Compute the Optimal Loop-Shaping Controller
[K,CL,GAM] = loopsyn(G,Gd);
GAM
K
%comparison
L = G*K;              % form the compensated loop L
figure

sigma(Gd,'b',L,'r--',{.1,100});
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times')
grid on
legend('Gd (target loop shape)','L (actual loop shape)');
%Analyzing the Shaped-Loop L, Closed-Loop T, and Sensitivity S
T = feedback(L,eye(1));
T.InputName = {'Force'};
S = eye(1)-T;

% SIGMA frequency response plots
figure

sigma(inv(S),'m',T,'g',L,'r--',Gd,'b',Gd/GAM,'b:',...
	Gd*GAM,'b:',{.1,100})
grid on
legend('1/\sigma(S) performance',...
	'\sigma(T) robustness',...
	'\sigma(L) open loop',...
	'\sigma(Gd) target loop shape',...
	'\sigma(Gd) \pm GAM(dB)');
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times')


set(findobj(gca,'Type','line','-not','Color','b'),'LineWidth',2);
h = findobj(gca,'Type','line','-not','Color','b');
set(h,'LineWidth',2);

% figure
% step(T,8)
% title('Responses to step commands');
% S = stepinfo(T)
figure
impulse(T,8)

title('Responses to impulse commands');
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times')