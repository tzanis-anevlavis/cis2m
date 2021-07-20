%% Simple system:
A = [1.0 0.2; ...
     0.0 1.0];
B = [0.22; 0.2];
S = Polyhedron('lb', [-5; -5], 'ub', [5; 5]);
U = Polyhedron('lb', -2, 'ub', 2);

%% Level 1:
disp("L=1");
cis1 = computeRCIS(A,B,S.A,S.b,U.A,U.b,[],[],[],0,1);
for L = 1:length(cis1)
    cis1(L).minHRep; 
    disp(cis1(L).volume);
end

%% Level 2:
disp("L=2");
cis2 = computeCIS(A,B,S.A,S.b,U.A,U.b,[],[],[],0,2);
for L = 1:length(cis2)
    cis2(L).minHRep; 
    disp(cis2(L).volume);
end

%% Level 3:
disp("L=3");
cis3 = computeCIS(A,B,S.A,S.b,U.A,U.b,[],[],[],0,3);
for L = 1:length(cis3)
    cis3(L).minHRep; 
    disp(cis3(L).volume);
end

%% Level 4:
disp("L=4");
cis4 = computeCIS(A,B,S.A,S.b,U.A,U.b,[],[],[],0,4);
for L = 1:length(cis4)
    cis4(L).minHRep; 
    disp(cis4(L).volume);
end

%% Level 5:
disp("L=5");
cis5 = computeCIS(A,B,S.A,S.b,U.A,U.b,[],[],[],0,5);
for L = 1:length(cis5)
    cis5(L).minHRep; 
    disp(cis5(L).volume);
end

%% Level 6:
disp("L=6");
cis6 = computeCIS(A,B,S.A,S.b,U.A,U.b,[],[],[],0,6);
for L = 1:length(cis6)
    cis6(L).minHRep; 
    disp(cis6(L).volume);
end

%% Plot:
close all;

figure;
hold on;
plot(S,'color','blue');
plot(cis6(1),'color','orange');
plot(cis5(1),'color','yellow');
plot(cis4(1),'color','green');
plot(cis3(1),'color','cyan');
plot(cis2(1),'color','gray');
plot(cis1,'color','white');
hold off;

%% %%%% Compute MCIS using MPT3 %%%%
% tic
% system = LTISystem('A',A,'B',B);
% mcis = system.invariantSet('X',S,'U', U, 'maxIterations',200);
% timeMCIS = toc;
% disp(mcis.volume);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%