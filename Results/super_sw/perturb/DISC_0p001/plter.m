
load w3_am4_n100.txt
load w3_am6_n100.txt
load w3_am8_n100.txt
load w3_upwind_n100.txt

plot(w3_upwind_n100(:,1),(w3_upwind_n100(:,4)-w3_upwind_n100(:,6))-(w3_upwind_n100_stationary(:,2)-w3_upwind_n100_stationary(:,6)),'linewidth',1.5)
hold on
plot(w3_upwind_n100(:,1),(w3_am8_n100(:,4)-w3_am8_n100(:,6))-(w3_am8_n100_stationary(:,2)-w3_am8_n100_stationary(:,6)),'linewidth',1.5)
plot(w3_am6_n100(:,1),(w3_am6_n100(:,4)-w3_am6_n100(:,6))-(w3_am6_n100_stationary(:,2)-w3_am6_n100_stationary(:,6)),'linewidth',1.5)
plot(w3_am4_n100(:,1),(w3_am4_n100(:,4)-w3_am4_n100(:,6))-(w3_am4_n100_stationary(:,2)-w3_am4_n100_stationary(:,6)),'linewidth',1.5)
grid on
xlabel('$x $','Fontsize',18,'Interpreter','latex')
ylabel('$h-h_{in} $','Fontsize',18,'Interpreter','latex')
legend('upind no wb','AM3','AM4','AM6')
set(legend,'FontSize',10);
str=['weno3_sup_DISCsmall_n100_error_h'];
savefig(str);
     s=struct('Version',1,'Format','pdf','Width','800','Height','200','Units','points','Resolution','300');
        %hgexport(gcf,str,s)
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,...
       'PaperPosition',[0 0 screenposition(3:4)],...
       'PaperSize',[screenposition(3:4)]);
        print('-dpdf',str)