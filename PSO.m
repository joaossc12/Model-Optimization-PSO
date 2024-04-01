function [Fbest, gbest, it] = PSO(costf, limites, solucaoInicial,parametros)
%OTIMIZAÇÃO POR ENXAME DE PARTICULAS ADAPTATIVO
%
%[fopt, xopt] = PSO(costf, limitesMax, parametros.NP)
%costf: funcao custo y = costf(x) a ser otimizada, onde x um vetor 1xD é a posicao de
%uma unica particula.
%limites: matrix 4xD que gere os limites do parametros do PSO, de modo que
%cada linha representaa posicao Maxima e minima e velocidades Maxima e
%minima respectimante.
%D: dimensao, ou quantidade de variaveis de desição.

it = 1;
tempoInicial = tic;





D = size(limites,2);

xMax= limites(1,:);
xMin = limites(2,:);

vMax = 0.2*(xMax-xMin);


% calcula posicao e velocidades iniciais
    for d = 1:D
        x(1:parametros.NP,d) = xMin(d) + rand(parametros.NP,1)*( xMax(d) - xMin(d) );
        v(1:parametros.NP,d) =  2*rand(parametros.NP,1)*vMax(d) - vMax(d);
    end
    if ~isempty(solucaoInicial)
        x(1:size(solucaoInicial,1),:) = max(min(solucaoInicial,xMax),xMin);
    end

% calcula o funcional custo
for n = 1:parametros.NP
    F(n) = costf(x(n, :));
end
tempoIt = toc(tempoInicial);

% determina Fbest e gbest
[C, I] = min(F);
Fbest = C;
gbest = x(I, :);

% pbest = x pois é a 1a iteracao
pbest = x;
Fpbest = F;

w = 0.9;

% atualiza a velocidade e posicao
G = gbest.*ones(parametros.NP,1);
v1 = rand(parametros.NP,D).*(parametros.cP*(pbest-x));
v2 = rand(parametros.NP,D).*(parametros.cS*(G-x));
v = w*v + v1 + v2;
v = min(vMax,v);
v = max(-vMax,v);
x = x + v;

x = min(xMax, x);
x = max(xMin, x);

while(1)
   %verifica criteiro de parada
    if ( it >= parametros.itMax )||( toc(tempoInicial)+tempoIt > parametros.tempoMax ) || (Fbest == 0)
        %infoOtm(it,Fbest,gbest,parametros.path); comentei
        fprintf('\t Tempo %9.4g',(toc(tempoInicial)));
        break;
    end
    
    if parametros.info>0
    %infoOtm(it,Fbest,gbest,parametros.path); comentei
    end
    
    % calcula o funcional custo
    for n = 1:parametros.NP
        F(n) = costf(x(n, :));
    end

    % determina Fbest e gbest
    [C, I] = min(F);
    if C<=Fbest
        Fbest = C;
        gbest = x(I, :);
    end
    
    
    i = F<Fpbest;
    pbest(i,:) = x(i,:);
    Fpbest(i) = F(i);
    
    
    %calculo da inercia (w) em funcao da distribuixao de distancias
    for n = 1:parametros.NP
        d(n) = 0;
        for j = 1:parametros.NP
            d(n) = d(n) + sqrt( sum( (x(n,:)-x(j,:)).^2, 'all' ));
        end
        d(n) = 1/parametros.NP*d(n);
    end
    
    dMax = max(d);
    dMin = min(d);
    
    f = (d(I)-dMin)/(dMax-dMin);
    w = 1/(1 + 1.5*exp(-2.6*f));
    
    % atualiza a velocidade e posicao
    G = gbest.*ones(parametros.NP,1);
    v1 = rand(parametros.NP,D).*(parametros.cP*(pbest-x));
    v2 = rand(parametros.NP,D).*(parametros.cS*(G-x));
    v = w*v + v1 + v2;
    v = min(vMax,v);
    v = max(-vMax,v);
    x = x + v;
    x = min(xMax, x);
    x = max(xMin, x);
    
    it = it+1;
end
    
end

