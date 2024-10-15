
load w3_am4_n100.txt
load w3_am8_n100.txt
load w3_upwind_n100.txt
load w3_am4_n100_st.txt
load w3_am8_n100_st.txt
load w3_upwind_n100_st.txt


%--eta
subplot(2,1,1)
plot(w3_upwind_n100(:,1),w3_upwind_n100(:,4)-w3_upwind_n100_st(:,4),'linewidth',1.5)
hold on
%plot(w3_upwind_n100(:,1),w3_am8_n100(:,4),'linewidth',1.5)
plot(w3_am8_n100(:,1),w3_am8_n100(:,4)-w3_am8_n100_st(:,4),'linewidth',1.5)
plot(w3_am4_n100(:,1),w3_am4_n100(:,4)-w3_am4_n100_st(:,4),'linewidth',1.5)
xlabel('$x $','Fontsize',18,'Interpreter','latex')
ylabel('$h-h_{in} $','Fontsize',18,'Interpreter','latex')
%ylabel('$\eta $','Fontsize',18,'Interpreter','latex')
grid on
legend('WENO3','GF-AM8','GF AM6')

subplot(2,1,2)
plot(w3_am8_n100(:,1),w3_am8_n100(:,4)-w3_am8_n100_st(:,4),'linewidth',1.5);
hold on
plot(w3_am4_n100(:,1),w3_am4_n100(:,4)-w3_am4_n100_st(:,4),'linewidth',1.5)

grid on
xlabel('$x $','Fontsize',18,'Interpreter','latex')
ylabel('$h-h_{in} $','Fontsize',18,'Interpreter','latex')
%ylabel('$\eta $','Fontsize',18,'Interpreter','latex')
legend('GF-AM8','GF-AM6')
set(legend,'FontSize',10);
 str=['weno3_sub_DISCsmall_n100_eta'];
 savefig(str);
      s=struct('Version',1,'Format','pdf','Width','800','Height','200','Units','points','Resolution','300');
         %hgexport(gcf,str,s)
         set(gcf,'Units','inches');
         screenposition = get(gcf,'Position');
         set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
         print('-dpdf',str)