# RadiobiologyModels
% Calculates NTCP for patients using dose volume histogram (DVH) data. 
% LKB equations can be found in Hanbook of Radiotherapy Physics: Theory &
% Practice Chpter 36 (36.3 - NTCP).

%clear all;close all;clc;
global inputfiles Pneu Frac DPF %(Pneu = pneumonitis, no pneumonitis. Frac = fraction vol irradited, DPF = dose/fraction)
inputfiles=dir('*.txt'); %DVH data for each patient in txt format
M = csvread('allboth.csv');
Pneu = M(:,29);
Frac = M(:,9);
DPF = M(:,4);

for i=1:length(inputfiles) %import patients
    fid(i)=fopen(inputfiles(i).name);
    inputfiles(i).values=textscan(fid(i),'%f%*f%f%[^\n\r]','delimiter',';','HeaderLines',33);
    fclose(fid(i));
end
n = 1.03;
TD50 = 24.5;
m = 0.37;

a0 = [n TD50 m]
options = optimset('MaxFunEvals', 100000,'MaxIter',10000,'Display','off');

a = fminsearch(@myfunc,a0, options) %fminsearch used for max likelihood estimation 
%N.B. does not seem to be working correctly


n = a(1);
TD50 = a(2);
m = a(3);

D = linspace(0,50,201); 
x = (D - TD50)/(m*TD50);

plot(D,normcdf(x));
hold on
plotPatients(n, TD50, m) %plots patients ntcp with fit parameters

n = 1.03;
TD50 = 24.5;
m = 0.37;

D = linspace(0,50,201);
x = (D - TD50)/(m*TD50);
plot(0, 0, 'Parent', AxesH);
plotPatients(n, TD50, m)%plots patients ntcp with literature parameters
% hold on

% plot(D,normcdf(x),'-b');

% legend('Moissenko et al (2003)','Fit','Location','SE')

xlabel('EUD (Gy)')
ylabel('NTCP')

function p = NTCP(n, TD50, m) % function to calculate NTCP
global inputfiles DPF 
p = zeros(length(inputfiles),1);
Deff = zeros(length(inputfiles),1);
f = Frac;
dpf = DPF;

for i=1:length(inputfiles)
    dose = cell2mat(inputfiles(i).values(:,1));
    vol = cell2mat(inputfiles(i).values(:,2));
    
    ab = 3.0;
    d = dpf(i);
    eqd2 = dose*((ab+d)/(ab+2)); %LQ 2Gy/fraction correction
   
    v = vol/sum(vol);
    dn = eqd2.^(1/n);
    frac = f(i).^n;
    
    
    Deff(i) = (sum(v.*dn))^n; %Effective Dose for each Patient

    t = (Deff(i) - (TD50/frac))/(m*(TD50/frac)) ;
    
    p(i) = normcdf(t); % NTCP calculation for each patient
end

end

function y = myfunc(v) % maxlike.lihood function - not working correctly
global inputfiles Pneu
ep = Pneu;

n = v(1);
TD50 = v(2);
m = v(3);

p = NTCP(n, TD50, m);
tmp = 0;
for i=1:length(inputfiles)
    tmp =  tmp + ep(i)*log(p(i)) + (1.0-ep(i))*log(1-p(i)); %max likelihood   
end
y = -tmp;
end

function plotPatients(n, TD50, m)%function calculates effective dose & ntcp 
global inputfiles Pneu DPF                   %again for patient plots
p = zeros(length(inputfiles),1);
Deff = zeros(length(inputfiles),1);

f = Frac;
dpf = DPF;
for i=1:length(inputfiles)
    dose = cell2mat(inputfiles(i).values(:,1));
    vol = cell2mat(inputfiles(i).values(:,2));
    
    ab = 3.0;
    d = dpf(i);
    
    eqd2 = dose*((ab+d)/(ab+2));
   
    v = vol/sum(vol);
    dn = eqd2.^(1/n);
    
    frac = f(i).^n;
    
    Deff(i) = (sum(v.*dn))^n;
    
    t = (Deff(i) - (TD50/frac))/(m*(TD50/frac));
    
    p(i) = normcdf(t);
end

numYes = 0; %Patients with Pneumonitis
numNo = 0; %Patients without Pneumonitis
for i=1:length(inputfiles)
    if Pneu(i) == 1
        numYes = numYes + 1;
    else
        numNo = numNo + 1;
    end
end
Dy = zeros(numYes,1); % patients with Pneu effective dose
Dn = zeros(numNo,1); % patients without Pneu effective dose

ntcpy = zeros(numYes,1); %ntcp with Pneu
ntcpn = zeros(numNo,1); %ntcp without Pneu

iy = 1;
in = 1;
for i=1:length(inputfiles)
    if Pneu(i) == 1
        Dy(iy) = Deff(i);
        ntcpy(iy) = p(i) ;
        iy = iy + 1;
    else
        Dn(in) = Deff(i);
        ntcpn(in) = p(i);
        in = in + 1;
    end
end
predicted = @(a, xdata) normcdf((xdata-a(1))/a(2))
a0 = [10, 0.5]
[ahat,r,J,cov,mse] = nlinfit(Deff,p,predicted,a0) % nonlinfit for nomcdf function (mu, sigma)
ci = nlparci(ahat,r,'Jacobian',J);
ci = nlparci(ahat,r,'covar',cov),

xplot = linspace(0,50,201);
yplot = normcdf((xplot-ahat(1))/ahat(2)); %NTCP fit 
plot(Dy, ntcpy, 'rx');
plot(Dn, ntcpn, 'ko');
%disp(ntcpy) 
%disp(ntcpn)
% disp(Dy)
% disp(Dn)
plot(xplot, yplot,'r'); %plot non-lin fit

str = {'x = Pneumonitis','o = No Pneumonitis'};
% str = {'x = Pneumonitis'};
text(10,0.6,str)
set(gca,'FontWeight','bold')
str1 = {'n = 0.8','m = 0.37', 'TD50 = 21.9'};
text(22,0.27,str1)
x = [19];
y = [0.45];
end
