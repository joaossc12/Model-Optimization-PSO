%Carrega dados 
%%
cd 'C:\Users\João Vitor\Desktop\UESC\IC\DadosSimul'
dados = csvread('bigdata.csv');
posicao = dados(:,1:3);
setpoint = dados(:,4:5);
U = setpoint';
tam = size(dados,1);
tempo =  0:0.5:((tam-1)*0.5);

X0 = zeros(7,1);
X = X0;
T = 0.5;
X1 = X0;

random = rand(1,8)*(1.3-0.7)+0.7;
H = [1, 1, 0.015, 10.913e-3, 1.550, 0.025, 0.02, 1].*random; %parametros
limites(1,:) = H*1.5;
limites(2,:) = H*0.5;

%% Pso MATLAB
tic()
options = optimoptions('particleswarm','SwarmSize',15,'MaxTime',10*60,'MaxStallIterations',50,'FunctionTolerance',0);

[xps, fxps,exitflag,output] = particleswarm(@(j)RMSE(j,posicao, tempo, U, X0, T),8,limites (2,:),limites (1,:),options);

TR2 = toc();

%% PSO Gabriel


tic()

parametros.itMax = inf;
parametros.tempoMax = 10*60;
parametros.NP = 15;
parametros.info = 1;
parametros.cP = 1.49; %PSO
parametros.cS = 1.49; %PSO
parametros.path = "";

[fopt, xopt,it] = PSO(@(h) RMSE(h,posicao, tempo, U, X1, T), limites, H, parametros);

TR1 =toc();

 
 %%
 %Simulação OTM

for i = 1:size(tempo, 2)
    [t, Xs] = ode45(@(t, X) ModeloDDMR(t, X, U(:, i), xopt), [tempo(i), tempo(i)+T], X0);
    X = Xs(end, :)'; % Atualiza as condições iniciais para a próxima iteração
    sim3(i, :) = X';
end

hold on
plot(posicao(:,1),posicao(:,2), 'r-');
plot(sim3(:,1),sim3(:,2), 'b-');
grid on;
xlabel('Eixo X');
ylabel('Eixo Y');
title('Comparação Modelo x Dados');
legend('Real','ModeloOTM');


