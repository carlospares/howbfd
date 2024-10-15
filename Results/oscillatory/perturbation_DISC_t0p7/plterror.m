load w5_am3_n150.txt
load w5_am4_n150.txt
load w5_am6_n150.txt
load w5_upwind_n150.txt


load w5_am3_n150_stationary.txt
load w5_am4_n150_stationary.txt
load w5_am6_n150_stationary.txt
load w5_upwind_n150_stationary.txt


% plot(w3_am3_n150(:,1),w3_am3_n150(:,3)-w3_am3_n150_stationary(:,2),'linewidth',1.5)
% hold on
% plot(w3_am4_n150(:,1),w3_am4_n150(:,3)-w3_am4_n150_stationary(:,2),'linewidth',1.5)
% plot(w3_am6_n150(:,1),w3_am6_n150(:,3)-w3_am6_n150_stationary(:,2),'linewidth',1.5)
% plot(w3_upwind_n150(:,1),w3_upwind_n150(:,3)-w3_upwind_n150_stationary(:,2),'linewidth',1.5)
% 

%------------weno3
% plot(w3_am3_n150(:,1),w3_am3_n150(:,3)-w3_am3_n150_stationary(:,3),'linewidth',1.5)
% hold on
% plot(w3_am4_n150(:,1),w3_am4_n150(:,3)-w3_am4_n150_stationary(:,3),'linewidth',1.5)
% plot(w3_am6_n150(:,1),w3_am6_n150(:,3)-w3_am6_n150_stationary(:,3),'linewidth',1.5)
% plot(w3_upwind_n150(:,1),w3_upwind_n150(:,3)-w3_upwind_n150_stationary(:,3),'linewidth',1.5)
% out=reference_p(:,3)-reference_p(:,2);
% out(500:1500)=0;
% %plot(reference_p(:,1),reference_p(:,3)-reference_stationary(:,3),'k','linewidth',1.5)
% plot(reference_p(:,1),out,'k','linewidth',1.5)

 %------------weno5

 plot(w5_am3_n150(:,1),w5_am3_n150(:,3)-w5_am3_n150_stationary(:,3),'linewidth',1.5)
 hold on
 plot(w5_am4_n150(:,1),w5_am4_n150(:,3)-w5_am4_n150_stationary(:,3),'linewidth',1.5)
 plot(w5_am6_n150(:,1),w5_am6_n150(:,3)-w5_am6_n150_stationary(:,3),'linewidth',1.5)
 plot(w5_upwind_n150(:,1),w5_upwind_n150(:,3)-w5_upwind_n150_stationary(:,3),'linewidth',1.5)
 out=reference_p(:,3)-reference_p(:,2);
 out(500:1500)=0;
 %plot(reference_p(:,1),reference_p(:,3)-reference_stationary(:,3),'k','linewidth',1.5)
 plot(reference_p(:,1),out,'k','linewidth',1.5)
 
grid on
xlabel('$x $','Fontsize',18,'Interpreter','latex')
ylabel('$u-u_{in} $','Fontsize',18,'Interpreter','latex')

legend('AM3','AM4','AM6','upind no wb','reference','Location','northwest')
set(legend,'FontSize',10);
%axis([-1 1 -1 3])

   str=['weno5_AM_DISC_n150'];
   savefig(str);
         s=struct('Version',1,'Format','pdf','Width','800','Height','200','Units','points','Resolution','300');
            %hgexport(gcf,str,s)
            set(gcf,'Units','inches');
            screenposition = get(gcf,'Position');
            set(gcf,...
           'PaperPosition',[0 0 screenposition(3:4)],...
           'PaperSize',[screenposition(3:4)]);
            print('-dpdf',str)