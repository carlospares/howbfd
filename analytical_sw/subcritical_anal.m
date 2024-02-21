clear h x r;

N=400;
dx=25/N;
x=[0.5/16:dx:25];
z = 0.05*sin(x-12.5).*exp(1-(x-12.5).^2);

%plot(x,z);

q0=24.; %(q(0,t)=q0)
h0=2.0;% (h(0,t)=h0=hL)
g=9.812;

a=ones(1,length(z));
b=z-q0^2./(2*g*h0^2) -h0;
c=zeros(1,length(z))+0.0;
d=ones(1,length(z))*q0^2/(2*g);

p=[a', b', c', d'];


for i=1:length(z)
    r=roots(p(i,:));
    h(i)=r(2);
    tmp1(i)=r(1);
    tmp2(i)=r(3);
end

sol=[ x', h', q0*ones(length(z),1)];
%sol=[h'];
%save('initial_sub_800.dat','sol','-ascii')
%fid = fopen('initial_sub_800.dat','w');
%fprintf(fid,'%2.16f \n',sol); %floating point with a width and precision of 30 numbers
%fclose(fid);