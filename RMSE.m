function Custo = RMSE(H,posicao, tempo, U, X, T)

for i = 1:size(tempo, 2)
    [t, Xs] = ode45(@(t, X) ModeloDDMR(t, X, U(:, i), H), [tempo(i), tempo(i)+T], X);
    X = Xs(end, :)'; % Atualiza as condições iniciais para a próxima iteração
    sim(i, :) = X';
end

simul = sim(:,1:3);
A = simul;
Erro = (posicao - simul).^2;
SumErro = sum(Erro, 'all');
EM = SumErro/size(posicao,1);

Custo = sqrt(EM);
disp(Custo);

end

