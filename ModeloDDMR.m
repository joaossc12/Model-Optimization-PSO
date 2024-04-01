function dX = ModeloDDMR(~,X,U, H)
%X -> Vetor de Estador [x, y, thetta, wd, we, Xed, Xee]
%U -> Entrada [ud; ue]
%H -> Parametros variáveis [Ic, Iw, b, Kce, Ra, Fs, Fk, alpha_s]

R = 0.033; % Raio das rodas do robô -OTP GABRIEL
L = 0.099; % Comprimento do semieixo das rodas do robô -OTP GABRIEL

Ic = H(1); % Momento de inércia total equivalente do robô -OTP GABRIEL
Iw = H(2); % Momento de inércia em torno do eixo de cada roda incluindo motor -OTP GABRIEL

b = H(3); % Coeficiente de atrito viscoso do motor CC -OTP GABRIEL
Kce = H(4);  % Constante de força contra-eletromotriz do motor CC -OTP GABRIEL
Kt = Kce; % Constante de torque do motor CC -OTP GABRIEL
Ra = H(5); % Resistencia de armadura total do robô
N = 19; % Relação de engrenagens do motor CC 

mc = 0.4; % Massa da plataforma do robô -OTP GABRIEL
mw = 0.2; % Massa da roda do robô -OTP GABRIEL

Fs = H(6); % Coistência de armadura do motor CC -OTP GABRIEL
Fk = H(7); % Diferença entre coeficiente de atrito estático e cinético -OTP GABRIEL 
alpha_s = H(8); % Constante de saturação do atrito com o solo 
alpha_k = alpha_s; % Constantes de saturação do atrito com o solo 

ki = 500; % Ganho Ki
kp = 100; % Ganho Kp


m = mc + 2*mw; %Massa do robô
It = Ic + 2*mw*L^2 + 2*Iw;
Ki = ki*eye(2,2);
Kp = kp*eye(2,2);
g = 9.81; %Aceleracao da gravidade


M_TRACO = [
    Iw+((R^2/(4*(L^2)))*(m*(L^2)+It)),      (R^2/(4*(L^2)))*(m*(L^2)-It);
    (R^2/(4*(L^2)))*(m*(L^2)-It),       Iw+((R^2/(4*(L^2)))*(m*(L^2)+It))];


B_TRACO = eye(2,2);

A_tal = (1/N)*[
    ((-Kce*Kt)/Ra)+b,           0;
    0,           ((-Kce*Kt)/Ra)+b];

B_tal = [
    Kt/Ra, 0;
    0, Kt/Ra];

A_tald = (M_TRACO\B_TRACO)*A_tal;

B_tald = (M_TRACO\B_TRACO)*B_tal;

Bx = [
    zeros(3,2);
    B_tald*Kp;
    eye(2,2)];

A_cin = (R/2)*[
    cos(X(3)),  cos(X(3));
    sin(X(3)),  sin(X(3));
    1/L,            -1/L];

Ax = [
    zeros(3,3),      A_cin,                      zeros(3,2);
    zeros(2,3),   (A_tald-B_tald*Kp),    (B_tald*Ki);
    zeros(2,3),     -eye(2,2)                   zeros(2,2)];

tal_a = [
    Fs*tanh(alpha_s*X(4))-Fk*tanh(alpha_k*X(4));
    Fs*tanh(alpha_s*X(5))-Fk*tanh(alpha_k*X(5))];


F_tald = -(M_TRACO\B_TRACO)*tal_a;

Fx = m*g*[
    zeros(3,1);
    F_tald;
    zeros(2,1)];

dX = Ax*X + Bx*U + Fx;



end

