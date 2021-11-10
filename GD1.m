clc
clear all
close all
%% GRADIEN DECENT TRAINING ALGORITHEM.
%   This programm develops as a course project for Fuzzy systems course by
%    Habibollah Naeimi.
disp(' GRADIEN DECENT TRAINING ALGORITHEM.');
disp('  This programm develops as a course project for Fuzzy systems course')
disp('   by Habibollah Naeimi.')
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% 1st Part: Parameter Setting.
%%  Parameters Initiating.
disp(' Parameters Initiating...');

M = 30;                                % M: Number of Rules.
%M = input('Please Enter your Rule Numbers(M):\n');
Alpha = 0.3;                           % Aplpha, The Training Ratio!
%Alpha = input('Please Enter Alpha:\n');
Q = 100;                               % Number of iteration for each point.
%PointR=input('Please Enter the Number of Iteration for each data point(q):\n');
epsilon = 0;                           % Final Error for Point Repeat Training
%e = input('Please Enter Desired Error for each point(epsilon!):\n');

InpuNum = 1;                           % Input Number.

y_Bar = zeros(M,1);                    % Preparing Initial Matrixes.
x_Bar = zeros(M,InpuNum);
Sigma = zeros(M,InpuNum);

disp(' Part 1: DONE!');
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% 2nd Part: Sampling.
%%  Calculating Samples.

Data_Pairs_Num = 300;                  % Number of Data Pairs.
%Data_Pairs_Num = input('Please Enter Number of Data Pairs:\n');
SAMPLES_Num = 500;                     % Number of Samples.
%SAMPLES_Num = input('Please Enter Number of Samples:\n');

y_aprx = zeros(1,SAMPLES_Num);         % Preparing Initial Matrixes.
SAMPLES = zeros(SAMPLES_Num,InpuNum+1);
iNSamples = 0.2:0.01:0.51;

for i=33:SAMPLES_Num+33+InpuNum
    iNSamples(i) = 0.2*iNSamples(i-31)/(1+(iNSamples(i-31)^10))+0.9*iNSamples(i-1);
end

iNSamples = iNSamples(33:end);
    
for i=1:SAMPLES_Num
    SAMPLES(i,:) = iNSamples(i:i+InpuNum);
end

Pairs = SAMPLES(1:Data_Pairs_Num,:);
disp(' Complete Sampling.');
%%  Primary Parameter Fixing Using Online initial Parameter Choosing.

x_Bar = Pairs(1:M,1:InpuNum);         
y_Bar = Pairs(1:M,end);
Sigma = repmat(((max(x_Bar)-min(x_Bar))/M),M,1);

disp(' Initial Parameters are Reasdy!');
disp(' ');
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% 3rd Part: Parameters Updating. 
%%  Updating.

z = zeros(M,1);                % Preparing Initial Matrixes.
iN_z = zeros(1,InpuNum);

 for p=1:size(Pairs,1)
   
     for l=1:M                 % Calculating z.
         for i=1:InpuNum
             iN_z(i) = exp(-(((Pairs(p,i)-x_Bar(l,i))/Sigma(l,i))^2));
         end
         z(l) = prod(iN_z);
     end
            
     b = sum(z);               % Calculating b.
     a = sum(y_Bar.*z);        % Calculating a.
     f = a/b;                  % Calculating f.

            for q=1:Q
                for l=1:M             % Updating Parameters.
                    y_Bar(l) = y_Bar(l)-Alpha*(f-Pairs(p,end))/b*z(l);
                    for i=1:InpuNum                    
                        x_Bar(l,i) = x_Bar(l,i)-Alpha*(f-Pairs(p,end))/b*(y_Bar(l)-f)*z(l)*(2*(Pairs(p,i)-x_Bar(l,i))/(Sigma(l,i)^2));
                        Sigma(l,i) = Sigma(l,i)-Alpha*(f-Pairs(p,end))/b*(y_Bar(l)-f)*z(l)*(2*((Pairs(p,i)-x_Bar(l,i))^2)/(Sigma(l,i)^3));
                    end
                end
                
                if (f-Pairs(p,end))<epsilon   % Error Checking.
                    break;
                end
            end
 end

disp(' Parameters Updated.');
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% 4th Part: Results.
%%  Results.
 
f = zeros(1,SAMPLES_Num);            % Preparing Initial Matrixes.
f(1:2) = SAMPLES(1:2,end);

for k=3:SAMPLES_Num
    
    for l=1:M                         % Calculating z.
        for i=1:InpuNum
            iN_z(i) = exp(-(((SAMPLES(k,i)-x_Bar(l,i))/Sigma(l,i))^2));                    
        end
        z(l) = prod(iN_z);
    end
           
    b = sum(z);                         % Calculating b.
    a = sum(y_Bar.*z);                  % Calculating a.
    f(k) = a/b;                         % Calculating f.
    y_aprx(k) = f(k);                   % Setting Approximation output.
             
end
%%  Plotting.

figure;
plot(SAMPLES(:,end));
hold on
plot(y_aprx,'r');
legend('Real Value','Predicted Value');
%%  Errors.

Error = SAMPLES(:,end)-y_aprx';
disp('Mean Square Error is:');
MSE = mse(Error)
disp('Mean Absolute Error is:');
MAE = mae(Error)






