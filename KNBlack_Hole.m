clear all;

Propotion=1.0;
Mass=Propotion*(1.99e+30)
Rs=2950.5*Propotion
RQ=0.4808/Propotion
amax=sqrt(0.25-RQ*RQ)
Energy=zeros(70,501);
AngukarMomentum=zeros(70,501);
Energy2=zeros(70,501);
AngukarMomentum2=zeros(70,501);
gamma=zeros(70,501);
gamma2=zeros(70,501);
mec2=(9.1e-31)*power(2.99792458e+8,2.0)
c=2.99792458e+8
const1=(9.388e+8)/(9.1e-31)
const3=amax/500
const7=0.001
const6=0.001

i_begin=1
for i=i_begin:70
    r=(0.9367+const7*(i-1));
    CoulombEnergy=(8.99e+9)*(1.60e-19)*(1.648824e+20)/(r*2950.5)
    for j=1:501
        a=const6*(j-1);
        const4=r*r-r+a*a+RQ*RQ;
        const8=r-2*a*a-2*RQ*RQ;
        const9=2*r*r-3*r+4*a*a+4*RQ*RQ;
          
        if const4>=0
           delta=sqrt(const4);
           Energy(i,j)=CoulombEnergy*((const8) - 2*a*delta)./(const9 - 4*a*delta);
           AngularMomentum(i,j)=(CoulombEnergy/c)*(a*(const8) - 2*(r*r-a*a)*delta)./(const9 - 4*a*delta);
           Energy2(i,j)=CoulombEnergy*((const8) + 2*a*delta)./(const9 + 4*a*delta);  %%順向
           AngularMomentum2(i,j)=(CoulombEnergy/c)*(a*(const8) + 2*(r*r-a*a)*delta)./(const9 + 4*a*delta);  %%順向
       end
    end
end
i_begin=i
gamma(:,:)=Energy(:,:)/mec2;
gamma2(:,:)=Energy2(:,:)/mec2;

figure(1)
for i=1:70;
     r=(0.9367+const7*(i-1));
     
     j=1:501;
     plot(const6*(j-1), Energy(i,j), 'r-');
     hold on
     plot(const6*(j-1), Energy2(i,j), 'g-');
     hold on
end
axis([0 0.03 -6e+9 0])
xlabel('Rotation a (Rs)')
ylabel('Energy E (J)');

figure(2)
for i=1:70  %% 33 33 電子群順轉時的負最大角動量出現
     r=(0.9367+const7*(i-1));
     
     j=1:501;
     plot(const6*(j-1), AngularMomentum(i,j), 'r-');
     hold on;
     plot(const6*(j-1), AngularMomentum2(i,j), 'g-');
     hold on;
end
axis([0 0.03 -40 5])
xlabel('Rotation a (Rs)');
ylabel('Angular Momentum L (Kg.m^2/s)');

for i=1:1000;
    alpha1(i)=0.001*i;  
    
    for j=1:1000;
        alpha2(j)=0.001*j;  

        RQ2(i,j)=(1-alpha1(i)-alpha2(j))*(1-alpha1(i)-alpha2(j))*RQ*RQ;
    end
end

for i1=1:70 
    r1(i1)=(0.9367+const7*(i1-1));
            
    r2(i1)=(0.9367+const7*(i1-1));
end 

alpha1_begin0=0.180
alpha2_begin0=0.280
alpha1_begin=round(alpha1_begin0/0.001)
alpha2_begin=round(alpha2_begin0/0.001)
i1_begin=1
n=0  %% 0.43
count=1
for i=alpha1_begin:1000;
    for j=alpha2_begin:1000;
        TotalMass=((Mass + (alpha1(i)+alpha2(j))*const1*(9.1e-31))*(2.99792458e+8));

        for i1=i1_begin:70 %% 33 33 電子群順轉時的負最大角動量出現
            for i2=1:70    %% 33 33 電子群順轉時的負最大角動量出現
                for j2=1:1  %% 250
                    const2=(alpha1(i)*AngularMomentum(i1,j2)*const1 + alpha2(j)*AngularMomentum2(i2,j2)*const1 + const6*(j-1)*(Mass)*(2.99792458e+8))/TotalMass;  %% 1836*
                    const8=1.0-4.0*const2^2-4.0*RQ2(i,j)+4.0*n*n;
                    const5=0.5*(1.0+sqrt(const8));
                    
                    if i1~=i2 && const8 >= 0 && const5 >= max(r1(i1),r2(i2)) %% abs(const2)<0.&& gamma(i,j)~=0
                       Solution_alpha1(count)=alpha1(i);
                       Solution_alpha2(count)=alpha2(j);
                       Solution_r1(count)=r1(i1);
                       Solution_r2(count)=r2(i2);
                       const5
                       AB(count)=const5;
                       count=count+1;
                   end
                end
            end
        end
    end
end
max(AB)

save Solution_alpha1.dat Solution_alpha1
save Soluttion_alpha2.dat Solution_alpha2
save Solution_r1.dat Solution_r1
save Solution_r2.dat Solution_r2
save AB.dat AB

figure(9)
i=1:count-1;
plot(Solution_r1(i), Solution_r2(i), 'r.');
hold on;
axis([0.935 1.00 0.935 1.00])
xlabel('Radius r1 (Rs)');
ylabel('Radius r2 (Rs)');






