%% Read txt file and chooose the size of your matrix
clear all
fid1 = fopen('output.txt');
line1=fscanf(fid1,'%s',1); line1=fscanf(fid1,'%s',1); line1=fscanf(fid1,'%s',1); line1=fscanf(fid1,'%s',1); line1=fscanf(fid1,'%s',1)
line1=fscanf(fid1,'%s',1)
line1=fscanf(fid1,'%s',1)

data1 = fscanf(fid1, '%e %e %e %e %e %e %e \n', [9 5000e3]); % inf]);
fclose(fid1);
disp(size(data1))
data2 = data1';
%% data2 = data(3000001:5000000,:);


Step = data2(:,1)'; % This is Step. 
Time = data2(:,2)'; % This is time; 
Temp = data2(:,3)'; % This is temperature;
Jx = data2(:,4)'; % This is Jx; 
Jy = data2(:,5)'; % This is Jy;
Jz = data2(:,6)'; % This is Jz;
pe = data2(:,7)';
ke = data2(:,8)';
disp('Size of Jx')
disp(size(Jx))
disp('Step');
disp(Step(1:15))

Time1 = Time - Time(1)*ones(1,length(Time));
Step1 = Time1;

%%---Input Information for different Geometries-----
%%--- Check these parameters carefully. ------------
T = 1000 % 1000;
disp('dt')
dt =  Step1(2) - Step1(1) % 0.550; % Si-bulk dt = 0.55 fs.   
s = 1; % print out Jx data every MD steps. 
a = 5.43071 % FOR Ar fcc. 5.4307; % lattice constance of Si. 
vol = (8*a)^3 %# *6*a*6*a; % 
%vol = 10082.5049187631 #9610.47; % 9962.1856; % 9781.3468;
taomax = 1000e3; %1.25e5; % 7e4; % 4e5; % 1000001;% 90910; % 3e5; % 1.5e5; %  % 90901 time_max = 500ps for bulk Si case.

%%---End: Input Information for different Geometries-----
% Time1 = Step.*dt; 
% Step1 = Time1; 
%----------If have Time already------------------
% Step2 = Step1; 
% % Time1 = Time - Time(1)*ones(length(Time),1);
% % Step1 = Time1;


%----------------------------
N = length(Jx);
J =[Jx Jy Jz];

kcal2j = 4184/6.022e23 ; %
fs2s = 1.0e-15; 
A2m = 1.0e-10; 
kB = 1.3806504e-23; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear G G_count Gx Gy Gz k kxx kyy kzz trapJJ trapJJx trapJJy trapJJz Rx

M = taomax;
    Gx = zeros(length(1:1:M),1);
    Gy = zeros(length(1:1:M),1);
    Gz = zeros(length(1:1:M),1);
    G = zeros(length(1:1:M),1);
%for m =1:1:M;% 

 %   Gx(m) = sum(Jx(m+1:N).*Jx(1:N-m))/(N-m);
 %   Gy(m) = sum(Jy(m+1:N).*Jy(1:N-m))/(N-m);
 %   Gz(m) = sum(Jz(m+1:N).*Jz(1:N-m))/(N-m);
 %   G(m) = (Gx(m)+Gy(m)+Gz(m))/3;
    %disp(m)
%end
%%%%-----------------------------------------------%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OCTAVE  Xcorr function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Rx lag] = xcorr(Jx, taomax, 'unbiased');
[Ry lag] = xcorr(Jy, taomax, 'unbiased');
[Rz lag] = xcorr(Jz, taomax, 'unbiased');

disp('xcorr finished')

Gx(1:taomax) = Rx(length(Rx)-taomax+1:length(Rx));
Gy(1:taomax) = Ry(length(Ry)-taomax+1:length(Ry));
Gz(1:taomax) = Rz(length(Rz)-taomax+1:length(Rz));
G  = (Gx+Gy+Gz)/3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dlmwrite("Gx.txt", Gx, "delimiter", "\\n");
%dlmwrite("Gy.txt", Gy, "delimiter", "\\n");
%dlmwrite("Gz.txt", Gz, "delimiter", "\\n");

%fid1 = fopen('Gx.txt');
%data1 = fscanf(fid1, '%e \n', [1 inf]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%------Integrate G, calculate k--------------------

trapJJ = zeros(length(1:taomax),1);
trapJJx = zeros(length(1:taomax),1);

trapJJy = zeros(length(1:taomax),1);

trapJJz = zeros(length(1:taomax),1);

trapJJ(1) = G(1);
trapJJx(1) = Gx(1);
trapJJy(1) = Gy(1);
trapJJz(1) = Gz(1);
for tao = 2:taomax; % 100001:1:200000;% 
    trapJJ(tao) = trapJJ(tao-1) +G(tao)*dt; % *tao ;
    trapJJx(tao) = trapJJx(tao-1) +Gx(tao)*dt; % *tao ;
    trapJJy(tao) = trapJJy(tao-1) +Gy(tao)*dt; % *tao ;
    trapJJz(tao) = trapJJz(tao-1) +Gz(tao)*dt; % *tao ;
end

disp('Finished trapJJ')

%%%%%%%%%%%%%%%%%%%%%For Heat Flux (divided by vol already in LAMMPS)%%%%
 
k = trapJJ*(kcal2j*kcal2j/fs2s/A2m)/kB/T/T*vol; % *s*dt;
kxx = trapJJx*(kcal2j*kcal2j/fs2s/A2m)/kB/T/T*vol; % *s*dt;
kyy = trapJJy*(kcal2j*kcal2j/fs2s/A2m)/kB/T/T*vol; % *s*dt;
kzz = trapJJz*(kcal2j*kcal2j/fs2s/A2m)/kB/T/T*vol; % *s*dt;

%%%%%%%%%%%%%%%%%%%%%For Heat Current (not divided by vol in LAMMPS)%%%%
%k = trapJJ*(kcal2j*kcal2j/fs2s/A2m)/kB/T/T/vol; % *s*dt; % /vol
%kxx = trapJJx*(kcal2j*kcal2j/fs2s/A2m)/kB/T/T/vol; % *s*dt;
%kyy = trapJJy*(kcal2j*kcal2j/fs2s/A2m)/kB/T/T/vol; % *s*dt;
%kzz = trapJJz*(kcal2j*kcal2j/fs2s/A2m)/kB/T/T/vol; %*s*dt;
disp('k finished')
csvwrite('k.csv',k)

%    save G.mat; % 
%    save Gx.mat 
%    save Gy.mat 
%    save Gz.mat
%  save k.mat
%  save kxx.mat
%  save kyy.mat
%  save kzz.mat
%%--------End of Calculation k----------------------

figure
plot(Step1(1:taomax)/1e3,G ,'LineWidth',4)
ylabel('HCACF');
xlabel('\tau (ps)')
legend(['Si-bulk-' num2str(T) 'K'])
print('HCACF-Si-bulk-test-','-dpng')
figure
plot(Step1(1:taomax)/1e3,Gx,'LineWidth',4)
ylabel('HCACF-x');
xlabel('\tau (ps)')
legend(['Si-bulk-' num2str(T) 'K'])
print('HCACF-x-Si-bulk-test','-dpng')

figure
plot(Step1(1:taomax)/1e3,Gy,'LineWidth',4)
ylabel('HCACF-y');
xlabel('\tau (ps)')
legend(['Si-bulk-' num2str(T) 'K'])
print('HCACF-y-Si-bulk-test','-dpng')

figure
plot(Step1(1:taomax)/1e3,Gz,'LineWidth',4)
ylabel('HCACF-z');
xlabel('\tau (ps)')
legend(['Si-bulk-' num2str(T) 'K'])
print('HCACF-z-Si-bulk-test','-dpng')

figure
plot(Step1(1:taomax)/1e3,G/G(1),'LineWidth',4)
ylabel('HCACF-normalized G/G(1)');
xlabel('\tau (ps)')
legend(['Si-bulk-' num2str(T) 'K'])
print('HCACF-normalize-Si-bulk','-dpng')

%%%%-----Plot k-----------------------

figure;
plot(Step1(1:taomax)/1e3, k,'LineWidth',4);
ylabel('k');
xlabel('\tau (ps)')
legend(['Si-bulk-' num2str(T) 'K'])
print('k-Si-bulk','-dpng');

figure;
plot(Step1(1:taomax)/1e3, kxx,'LineWidth',4);
ylabel('k-x');
xlabel('\tau (ps)')
legend(['Si-bulk-' num2str(T) 'K-output.txt'])
print('k-x-Si-bulk','-dpng');

figure;
plot(Step1(1:taomax)/1e3, kyy,'LineWidth',4);
ylabel('k-y');
xlabel('\tau (ps)')
legend(['Si-bulk-' num2str(T) 'K'])
print('k-y-Si-bulk','-dpng');

figure;
plot(Step1(1:taomax)/1e3, kzz,'LineWidth',4);
ylabel('k-z');
xlabel('\tau (ps)')
legend(['Si-bulk-' num2str(T) 'K'])
print('k-z-Si-bulk','-dpng');

