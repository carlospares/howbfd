%cd weno7_upwind
cd weno3_upwind/
load out25;
load out50;
load out100;
load out200;
load out400;
%load out800;

N=[25 50 100 200 400];% 800];
dx=25./N;

error1=[norm(out25(:,2)-out25(:,4),1)*dx(1), norm(out50(:,2)-out50(:,4),1)*dx(2), norm(out100(:,2)-out100(:,4),1)*dx(3), norm(out200(:,2)-out200(:,4),1)*dx(4), norm(out400(:,2)-out400(:,4),1)*dx(5)];%, norm(out800(:,2)-out800(:,4),1)*dx(6)];
%error2=[norm(out25(:,2)-out25(:,4),2)*sqrt(dx(1)), norm(out50(:,2)-out50(:,4),2)*sqrt(dx(2)), norm(out100(:,2)-out100(:,4),2)*sqrt(dx(3)), norm(out200(:,2)-out200(:,4),2)*sqrt(dx(4)), norm(out400(:,2)-out400(:,4),2)*sqrt(dx(5)), norm(out800(:,2)-out800(:,4),2)*sqrt(dx(6))];

%sum2=[sum(out25(:,2).^2), sum(out50(:,2).^2),sum(out100(:,2).^2),sum(out200(:,2).^2),sum(out400(:,2).^2) ,sum(out800(:,2).^2)];
%error2_rel=error2./sqrt(sum2);

rate1=log2(error1(2:end)./error1(1:end-1))./log2(dx(2:end)./dx(1:end-1))
%rate2=log2(error2(2:end)./error2(1:end-1))./log2(dx(2:end)./dx(1:end-1))
%rate2_rel=log2(error2_rel(2:end)./error2_rel(1:end-1))./log2(dx(2:end)./dx(1:end-1))

%error1_inout_bc7=error1;
%error2_inout_bc7=error2;
