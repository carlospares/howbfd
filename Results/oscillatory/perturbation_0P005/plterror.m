load w5_am3_n250.txt
load w5_am4_n250.txt
load w5_am6_n250.txt
load w5_upwind_n250.txt


load w5_am3_n250_stationary.txt
load w5_am4_n250_stationary.txt
load w5_am6_n250_stationary.txt
load w5_upwind_n250_stationary.txt


% plot(w3_am3_n150(:,1),w3_am3_n150(:,3)-w3_am3_n150_stationary(:,2),'linewidth',1.5)
% hold on
% plot(w3_am4_n150(:,1),w3_am4_n150(:,3)-w3_am4_n150_stationary(:,2),'linewidth',1.5)
% plot(w3_am6_n150(:,1),w3_am6_n150(:,3)-w3_am6_n150_stationary(:,2),'linewidth',1.5)
% plot(w3_upwind_n150(:,1),w3_upwind_n150(:,3)-w3_upwind_n150_stationary(:,2),'linewidth',1.5)
% 

%------------weno3
% plot(w3_am3_n250(:,1),w3_am3_n250(:,3)-w3_am3_n250_stationary(:,3),'linewidth',1.5)
% hold on
% plot(w3_am4_n250(:,1),w3_am4_n250(:,3)-w3_am4_n250_stationary(:,3),'linewidth',1.5)
% plot(w3_am6_n250(:,1),w3_am6_n250(:,3)-w3_am6_n250_stationary(:,3),'linewidth',1.5)
% plot(w3_upwind_n250(:,1),w3_upwind_n250(:,3)-w3_upwind_n250_stationary(:,3),'linewidth',1.5)
% %out=reference_p(:,3)-reference_p(:,2);
% %out(500:1500)=0;
% %%plot(reference_p(:,1),reference_p(:,3)-reference_stationary(:,3),'k','linewidth',1.5)
% %plot(reference_p(:,1),out,'k','linewidth',1.5)

 %------------weno5

 plot(w5_am3_n250(:,1),w5_am3_n250(:,3)-w5_am3_n250_stationary(:,3),'linewidth',1.5)
 hold on
 plot(w5_am4_n250(:,1),w5_am4_n250(:,3)-w5_am4_n250_stationary(:,3),'linewidth',1.5)
 plot(w5_am6_n250(:,1),w5_am6_n250(:,3)-w5_am6_n250_stationary(:,3),'linewidth',1.5)
 plot(w5_upwind_n250(:,1),w5_upwind_n250(:,3)-w5_upwind_n250_stationary(:,3),'linewidth',1.5)
%  out=reference_p(:,3)-reference_p(:,2);
%  out(500:1500)=0;
%  %plot(reference_p(:,1),reference_p(:,3)-reference_stationary(:,3),'k','linewidth',1.5)
%  plot(reference_p(:,1),out,'k','linewidth',1.5)
 
grid on
xlabel('$x $','Fontsize',18,'Interpreter','latex')
ylabel('$u-u_{in} $','Fontsize',18,'Interpreter','latex')

%legend('AM3','AM4','AM6','upind no wb','reference','Location','northwest')
legend('WENO5 GF-AM3','WENO5 GF-AM4','WENO5 GF-AM6','WENO5','Location','northwest')
set(legend,'FontSize',10);
axis([-1 1 0 20*10^-3])

     str=['weno5_AM_error_pert0p005'];
   