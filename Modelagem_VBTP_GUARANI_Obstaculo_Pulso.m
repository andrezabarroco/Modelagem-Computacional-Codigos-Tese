
 %===========================INSTITUTO MILITAR DE ENGENHARIA RJ===============================================
 %Autora: Andreza Barroco de Jesus
%Modelo computacional da VBTP GUARANI trafegando por obstáculo do tipo
%pulso, conforme resolução número 600 do CONTRAN
%Realiza análise da dinâmica vertical e do conforto.
%Análise através de índices de conforto como aceleração RMS e VDV.
%Otimização por Enxame de Partículas PSO para os parâmetros de rigidez e
%amortecimento da suspensão do tipo MacPherson.
 %==========================================================================

close all;clear all;clc
disp('Modelo computacional da VBTP GUARANI atravessando obstáculo conforme a Norma n 600 CONTRAN, realizando otimização através do algoritmo de PSO')

%% 1
%% ------------------MODELO DE PERFIL DO OBSTÁCULO----------------------
L = 350; % Comprimento total da estrada
Bd = 0.1; % Passo da discretização
% Intervalo de percurso observado (m)
x = 0:Bd:L;
 % Velocidade de percurso (m/s)
 v=20;
 % Tempo de percurso observado (s)
 t=x/v;
%% ------------------ DEFINIÇÃO DO PERFIL DO QUEBRA-MOLAS ----------------------
 % Formula da Entrada Lombada (m)--> Conforme Resolução no600 do CONTRAN
  h=[zeros(1,100) 0.1*sin(pi*(v*t(101:138)-x(101))/(x(138)-x(101))).^2 zeros(1,length(t)-138)];

%% ------------------- CÁLCULO DA VELOCIDADE E ACELERAÇÃO DA ESTRADA -------------------
v_estrada = gradient(h, Bd); % Velocidade da estrada (m/s)
a_estrada = gradient(v_estrada, Bd); % Aceleração da estrada (m/s²)

%% ----------------- DADOS DO MODELO DO VEÍCULO ----------------------
ms = 14247.92 / 2; % Massa suspensa (kg)
mp = 415; % Massa não suspensa (kg)
kp = 1250000; % Rigidez do pneu (N/m)

% Parâmetros de otimização
k1 = 395000; % Rigidez suspensão dianteira
k2 = 550000; % Rigidez suspensão intermediária
k3 = 458000; % Rigidez suspensão traseira

br = 27692.4; % Coeficiente de amortecimento Rebound
bB = 10000; % Rigidez amortecimento Bump

a1 = 1.7746; % eixo dianteiro
a2 = 0.0746; % eixo intermediário 
a3 = 1.9254; % eixo traseiro

Iyy = 52729.5; %Inércia eixo y

% ------------- Resultado das respostas dinâmicas (Saídas) -----------------
                % Montar a matriz de estado
        A1 = [0 1 0 0; -(kp + k1) / mp -br / mp k1 / mp br / mp; 0 0 0 1; k1 / ms br / ms -k1 / ms -br / ms];
        A2 = [0 1 0 0; -(kp + k2) / mp -bB / mp k2 / mp bB / mp; 0 0 0 1; k2 / ms bB / ms -k2 / ms -bB / ms];
        A3 = [0 1 0 0; -(kp + k3) / mp -br / mp k3 / mp br / mp; 0 0 0 1; k3 / ms br / ms -k3 / ms -br / ms];
        A = A1 + A2 + A3;

        C1 = [k1 / ms br / ms -k1 / ms -br / ms; 1 0 -1 0; -kp 0 0 0];
        C2 = [k2 / ms bB / ms -k2 / ms -bB / ms; 1 0 -1 0; -kp 0 0 0];
        C3 = [k3 / ms br / ms -k3 / ms -br / ms; 1 0 -1 0; -kp 0 0 0];
        C = C1 + C2 + C3;

        B = [0; kp / mp; 0; 0];
        D = [0; 0; kp];

        %% ---------- Cálculo da solução com a Entrada Perfil de Estrada ------------
        sys1 = ss(A1, B, C1, D);
        y1 = lsim(sys1, h, t);
        sys2 = ss(A2, B, C2, D);
        y2 = lsim(sys2, h, t);
        sys3 = ss(A3, B, C3, D);
        y3 = lsim(sys3, h, t);

        sys = ss(A, B, C, D);
        y = lsim(sys, h, t);

        %% ------------- Resultado das respostas dinâmicas (Saídas) -----------------
        % Aceleração vertical do chassi (m/s)
        As = y(:, 1);
        aceleracao_rms = sqrt(mean(As.^2)); % Aceleração RMS
        a_rms_before = aceleracao_rms;
        %% ------------- Resultado das respostas dinâmicas (Saídas) -----------------
        % Deflexão da suspensão (m)
        Ds1 = y1(:, 2);
        Ds2 = y2(:, 2);
        Ds3 = y3(:, 2);
        
        % Força dinâmica da roda (N)
        F1 = y1(:, 3);
        F2 = y2(:, 3);
        F3 = y3(:, 3);


        disp ('Ds1,Ds2,Ds3')
        disp (max(Ds1))
        disp (max(Ds2))
        disp (max(Ds3))
        disp ('F1,F2,F3')
        disp (max(F1))
        disp (max(F2))
        disp (max(F3))
        disp('Aceleração')
        disp(max(As))

%% Gráficos de Força Dinâmica da Roda
figure(11) 
hold on;
plot(t, F1, 'DisplayName', 'F1', 'LineWidth', 1.5);
plot(t, F2, 'DisplayName', 'F2', 'LineWidth', 1.5);
plot(t, F3, 'DisplayName', 'F3', 'LineWidth', 1.5);
title('Força Dinâmica das Rodas ao Longo do Tempo');
xlabel('t (s)');
ylabel('F(N)');
legend('show');
grid on;

%% ------------------- GERAÇÃO DE GRÁFICOS -------------------------
%% Gráfico 3: Deslocamento Vertical do Chassi ao Longo do Tempo
Z_chassi = cumtrapz(t, As);  
figure(3)
plot(t, Z_chassi, 'LineWidth', 1.5);
title('Deslocamento Vertical do Chassi ao Longo do Tempo', 'FontSize', 14);
xlabel('t (s)', 'FontSize', 12);
ylabel('Z (m)', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12); % Aumenta o tamanho da fonte dos eixos

disp('Bounce')
disp(max(Z_chassi))
%% Gráfico 1: Perfil da Estrada
figure(1) 
plot(x, h, 'LineWidth', 1.5);
title('Perfil da Estrada','FontSize', 14);
xlabel('x (m)','FontSize', 12);
ylabel('h (m)','FontSize', 12);
grid on;
set(gca, 'FontSize', 12); % Aumenta o tamanho da fonte dos eixos

%% Gráfico 2: Aceleração do Chassi
figure(2) 
plot(t, As, 'LineWidth', 1.5);
title('Aceleração Vertical do Chassi ao Longo do Tempo','FontSize', 14);
xlabel('t (s)','FontSize', 12);
ylabel('a (m/s²)','FontSize', 12);
grid on;
set(gca, 'FontSize', 12); % Aumenta o tamanho da fonte dos eixos


%% Gráficos de Deflexão da Suspensão
figure(4) 
hold on;
plot(t, Ds1, 'DisplayName', 'Ds1', 'LineWidth', 1.5);
plot(t, Ds2, 'DisplayName', 'Ds2', 'LineWidth', 1.5);
plot(t, Ds3, 'DisplayName', 'Ds3', 'LineWidth', 1.5);
title('Deflexão da Suspensão ao Longo do Tempo','FontSize', 14);
xlabel('t (s)');
ylabel('Ds (m)');
legend('show');
grid on;
set(gca, 'FontSize', 12); % Aumenta o tamanho da fonte dos eixos

%% Gráfico do perfil da estrada
figure(7)
subplot(2,1,1);
plot(x, h);
title('Perfil da Estrada','FontSize', 14);
xlabel('x (m)','FontSize', 12);
ylabel('h (m)','FontSize', 12);
grid on;


%% Gráfico da aceleração da estrada
subplot(2,1,2);
plot(x, a_estrada);
title('Aceleração da Estrada','FontSize', 14);
xlabel('x (m)','FontSize', 12);
ylabel('a (m/s²)','FontSize', 12);
grid on;
% Ajuste do layout
sgtitle('Análise do Perfil e Aceleração da Estrada'); % Título geral
set(gca, 'FontSize', 12); % Aumenta o tamanho da fonte dos eixos

%% Pitch
% Calculando o momento de pitch com base nas forças nas suspensões
M_pitch = F1 * a1 + F2 * a2 + F3 * a3;  % Momentos devido às forças nas suspensões

% Calculando o ângulo de pitch em radianos
pitch_angle = M_pitch / Iyy;

% Plotar o ângulo de pitch
figure(123)
plot(t, pitch_angle);
xlabel('t (s)');
ylabel('θ (rad)','FontSize', 14);
title('Movimento de Pitch');
grid on;
set(gca, 'FontSize', 12); % Aumenta o tamanho da fonte dos eixos


%% Frequência
% Perfil do obstáculo
N = length(h); % Número de pontos no perfil da estrada
dx = L / N; % Intervalo entre os pontos 

% Aplicar a Transformada de Fourier ao perfil da estrada
H = fft(h); % FFT do perfil da estrada

% Frequências associadas ao perfil da estrada
frequencias = (0:N-1) / (N * dx); % Frequências espaciais

% Magnitude da Transformada de Fourier 
magnitudes = abs(H);

% Plotar a magnitude da transformada de Fourier
figure;
plot(frequencias, magnitudes);
title('Análise de Frequência do Perfil da Estrada');
xlabel('Frequência (ciclos por metro)');
ylabel('Magnitude');

% Frequência máxima no perfil da estrada (em ciclos por metro)
f_max_estrada = max(frequencias);

% Converter para Hz
f_max_estrada_Hz = (f_max_estrada * v) / L;

disp(['Frequência máxima de excitação (Hz): ', num2str(f_max_estrada_Hz)]);

% Frequência natural para cada suspensão
f1 = (1/(2*pi)) * sqrt(k1/ms); % Frequência da suspensão dianteira (Hz)
f2 = (1/(2*pi)) * sqrt(k2/ms); % Frequência da suspensão intermediária (Hz)
f3 = (1/(2*pi)) * sqrt(k3/ms); % Frequência da suspensão traseira (Hz)

disp(['Frequência natural da suspensão dianteira (Hz): ', num2str(f1)]);
disp(['Frequência natural da suspensão intermediária (Hz): ', num2str(f2)]);
disp(['Frequência natural da suspensão traseira (Hz): ', num2str(f3)]);

%% Parâmetros do assento
m_assento = 35; % Massa do assento (kg)
k_assento = 15000; % Rigidez do assento (N/m)
b_assento = 10000 ; % Amortecimento do assento (N.s/m)

% Parâmetros dos Segmentos
M_1 = 4.17;   % Massa da Cabeça e o Pescoço (kg)
M_2 = 15;     % Massa do Torso Superior (kg)
M_3 = 5.5;    % Massa do Torso Inferior (kg)
M_4 = 36;     % Massa das Coxas e Pélvis (kg)

B_1 = 310;    % Amortecimento do Pescoço (N.s/m)
B_2 = 909.1 + 200;  % Amortecimento do Torso Superior (N.s/m)
B_3 = 330;    % Amortecimento do Torso Inferior (N.s/m)
B_4 = 2475;   % Amortecimento das Coxas e Pélvis (N.s/m)

K_1 = 16699;  % Rigidez do Pescoço (N/m)
K_2 = 10000;  % Rigidez do Torso Superior (N/m)
K_3 = 20000;  % Rigidez do Torso Inferior (N/m)
K_4 = 49340;  % Rigidez das Coxas e Pélvis (N/m)

%% Matriz de estado para o sistema massa-mola-amortecedor (assento)
A_assento = [0 1; -k_assento/m_assento -b_assento/m_assento];
B_assento = [0; 1];
C_assento = [1 0]; % Saída: deslocamento do assento
D_assento = 0;
a_chassi = As;
% Criar sistema de espaço de estados
sys_assento = ss(A_assento, B_assento, C_assento, D_assento);


% Resolver a resposta dinâmica do assento usando a aceleração do chassi como entrada
Z_assento = lsim(sys_assento, m_assento * a_chassi, t);
t_span=t';

% Matrizes de Estado para Cada Segmento
A_1 = [0, 1; -K_1/M_1, -B_1/M_1];  % Matriz para a cabeça
A_2 = [0, 1; -K_2/M_2, -B_2/M_2];  % Matriz para o torso superior
A_3 = [0, 1; -K_3/M_3, -B_3/M_3];  % Matriz para o torso inferior
A_4 = [0, 1; -K_4/M_4, -B_4/M_4];  % Matriz para as coxas/pélvis

B_f_1 = [0; 1/M_1];  % Força aplicada no segmento 1 (cabeça)
B_f_2 = [0; 1/M_2];  % Força aplicada no segmento 2 (torso superior)
B_f_3 = [0; 1/M_3];  % Força aplicada no segmento 3 (torso inferior)
B_f_4 = [0; 1/M_4];  % Força aplicada no segmento 4 (coxa/pélvis)

% Definindo a matriz
C = [1, 0]; 
D = 0;      

% Sistema de Estado para Cada Segmento (utilizando o comando ss)
sys_1 = ss(A_1, B_f_1, C, D);  % Sistema para o segmento 1 (cabeça)
sys_2 = ss(A_2, B_f_2, C, D);  % Sistema para o segmento 2 (torso superior)
sys_3 = ss(A_3, B_f_3, C, D);  % Sistema para o segmento 3 (torso inferior)
sys_4 = ss(A_4, B_f_4, C, D);  % Sistema para o segmento 4 (coxa/pélvis)
x1 = lsim(sys_1, Z_assento, t_span);  % Resposta para o segmento 1 (cabeça)
x2 = lsim(sys_2, Z_assento, t_span);  % Resposta para o segmento 2 (torso superior)
x3 = lsim(sys_3, Z_assento, t_span);  % Resposta para o segmento 3 (torso inferior)
x4 = lsim(sys_4, Z_assento, t_span);  % Resposta para o segmento 4  (coxa/pélvis)

% Extrair a deflexão e a aceleração 
deflexao1 = x1(:, 1);  % Deflexão da cabeça
deflexao2 = x2(:, 1);  % Deflexão do torso superior
deflexao3 = x3(:, 1);  % Deflexão do torso inferior
deflexao4 = x4(:, 1);  % Deflexão das coxas/pélvis

figure
plot(t, deflexao1, t, deflexao2, t, deflexao3, t, deflexao4);
title('Deflexão dos Segmentos do Corpo');
legend('Cabeça', 'Torso Superior', 'Torso Inferior', 'Coxas/Pélvis');
grid on
xlabel('t (s)');
ylabel('D (m)');

%% Encontrar o máximo deslocamento de cada parte
max_deslocamento_cabeca = max(deflexao1);
max_deslocamento_torso_superior = max(deflexao2);
max_deslocamento_torso_inferior = max(deflexao3);
max_deslocamento_pernas = max(deflexao4);

% Exibir os máximos deslocamentos
disp(['Máximo Deslocamento da Cabeça: ', num2str(max_deslocamento_cabeca), ' m']);
disp(['Máximo Deslocamento do Torso Superior: ', num2str(max_deslocamento_torso_superior), ' m']);
disp(['Máximo Deslocamento do Torso Inferior: ', num2str(max_deslocamento_torso_inferior), ' m']);
disp(['Máximo Deslocamento das Pernas: ', num2str(max_deslocamento_pernas), ' m']);

%% PSO PARAMETERS
n_particles = 30; % Número de partículas
n_iterations = 100; % Número de iterações
max_iter = n_iterations;
lower_bound = [200000, 200000, 200000, 10000, 10000]; % Limites inferiores
upper_bound = [800000, 800000, 800000, 50000, 50000]; % Limites superiores
dim = length(lower_bound); % Dimensões do espaço de busca (k1, k2, k3, br, bB)

%% Antes da otimização
% Inicialização das partículas
positions = lower_bound + (upper_bound - lower_bound) .* rand(n_particles, dim);
velocities = zeros(n_particles, dim);
personal_best_positions = positions;
personal_best_scores = inf(n_particles, 1);
global_best_score = inf;
global_best_position = zeros(1, dim);

% Parâmetros definidos
w1 = 0.5;  
w2 = 0.5;  
C1 = 2; % 1 a 2
C2 = 2;

%% Loop de iterações do PSO
for iter = 1:n_iterations
    % Cálculo do Fitness Final para cada partícula
    for i = 1:n_particles
max_distance = norm([500000, 500000, 500000, 25000, 25000] - upper_bound);

% Cálculo do fitness inicial
fitness_final_a = sum((positions(i, :) - [500000, 500000, 500000, 25000, 25000]).^2) / max_distance;
fitness_final = w1*fitness_final_a+w2*a_rms_before;

        % Avaliação da função de fitness da partícula
        if fitness_final < personal_best_scores(i)
            personal_best_scores(i) = fitness_final;
            personal_best_positions(i, :) = positions(i, :);
        end

        if fitness_final < global_best_score
            global_best_score = fitness_final;
            global_best_position = positions(i, :);
        end
    end

    % Atualização das velocidades e posições
    for i = 1:n_particles
        r1 = rand(1, dim);
        r2 = rand(1, dim);
        velocities(i, :) = 0.5 * velocities(i, :) + ...
            C1 * r1 .* (personal_best_positions(i, :) - positions(i, :)) + ...
            C2 * r2 .* (global_best_position - positions(i, :));
        positions(i, :) = positions(i, :) + velocities(i, :);

        % Verificação dos limites
        positions(i, :) = max(min(positions(i, :), upper_bound), lower_bound);
    end

    %% Definindo as matrizes K1, K2 e Z
    k1_range = unique(positions(:, 1)); 
    k2_range = unique(positions(:, 2)); 
    [K1, K2] = meshgrid(k1_range, k2_range); 
    Z = griddata(positions(:, 1), positions(:, 2), positions(:, 3), K1, K2); 

    % Gráfico de Superfície
    figure(81);
    surf(K1, K2, Z); 
    title('Superfície das Posições das Partículas', 'FontSize', 14);
    xlabel('k1', 'FontSize', 12); ylabel('k2', 'FontSize', 12); zlabel('k3', 'FontSize', 12);
    shading interp; 
    colorbar; 
    drawnow;
end

% Exibição dos resultados antes da otimização
disp('Comparativo dos resultados antes da otimização:');
disp(['Melhor score antes: ', num2str(global_best_score)]);
disp(['Melhor posição antes: ', num2str(global_best_position)]);
disp(['Aceleração RMS antes: ', num2str(a_rms_before)]);

% Valores otimizados da suspensão
k1_opt = global_best_position(1);
k2_opt = global_best_position(2);
k3_opt = global_best_position(3);
br_opt = global_best_position(4);
bB_opt = global_best_position(5);

A1_opt = [0 1 0 0; -(kp+k1_opt)/mp -br_opt/mp k1_opt/mp br_opt/mp; 0 0 0 1; k1_opt/ms br_opt/ms -k1_opt/ms -br_opt/ms];
A2_opt = [0 1 0 0; -(kp+k2_opt)/mp -bB_opt/mp k2_opt/mp bB_opt/mp; 0 0 0 1; k2_opt/ms bB_opt/ms -k2_opt/ms -bB_opt/ms];
A3_opt = [0 1 0 0; -(kp+k3_opt)/mp -br_opt/mp k3_opt/mp br_opt/mp; 0 0 0 1; k3_opt/ms br_opt/ms -k3_opt/ms -br_opt/ms];
A_opt = A1_opt + A2_opt + A3_opt;

sys_opt = ss(A_opt, [0; kp/mp; 0; 0], [k1_opt/ms br_opt/ms -k1_opt/ms -br_opt/ms], 0);
y_opt = lsim(sys_opt, h, t);

As_opt = y_opt(:, 1); % Aceleração vertical após otimização
a_rms_after = rms(As_opt); % Cálculo da aceleração RMS após otimização

        disp('Aceleração Otimizada')
        disp(max(As_opt))

% Frequência natural para cada suspensão
f1_opt = (1/(2*pi)) * sqrt(k1_opt/ms); % Frequência da suspensão dianteira (Hz)
f2_opt = (1/(2*pi)) * sqrt(k2_opt/ms); % Frequência da suspensão intermediária (Hz)
f3_opt = (1/(2*pi)) * sqrt(k3_opt/ms); % Frequência da suspensão traseira (Hz)

disp(['Frequência natural da suspensão dianteira (Hz): ', num2str(f1_opt)]);
disp(['Frequência natural da suspensão intermediária (Hz): ', num2str(f2_opt)]);
disp(['Frequência natural da suspensão traseira (Hz): ', num2str(f3_opt)]);

% Parâmetros definidos
w1 = 0.5;  
w2 = 0.5;  
C1 = 2; % 1 a 2
C2 = 2;  

% Fitness Final 
fitness_final_opt = w1 * fitness_final + w2 * a_rms_after;   

% Inicialização das partículas
positions_opt = lower_bound + (upper_bound - lower_bound) .* rand(n_particles, dim);
velocities_opt = zeros(n_particles, dim);
personal_best_positions_opt = positions;
personal_best_scores_opt = inf(n_particles, 1);
global_best_score_opt = inf;
global_best_position_opt = zeros(1, dim);

for i = 1:n_particles
    % Avaliação da função de fitness da partícula
    if fitness_final_opt < personal_best_scores_opt(i)
        personal_best_scores_opt(i) = fitness_final_opt;
        personal_best_positions_opt(i, :) = positions_opt(i, :);
    end

    if fitness_final_opt < global_best_score_opt
        global_best_score_opt = fitness_final_opt;
        global_best_position_opt = positions_opt(i, :);
    end
end

% Atualização das velocidades e posições
for i = 1:n_particles
    r1 = rand(1, dim);
    r2 = rand(1, dim);
    velocities_opt(i, :) = 0.5 * velocities_opt(i, :) + ...
        C1 * r1 .* (personal_best_positions_opt(i, :) - positions_opt(i, :)) + ...
        C2 * r2 .* (global_best_position_opt - positions_opt(i, :));
    positions_opt(i, :) = positions_opt(i, :) + velocities_opt(i, :);

    % Verificação dos limites
    positions_opt(i, :) = max(min(positions_opt(i, :), upper_bound), lower_bound);

    %% Visualização das partículas ao longo do tempo
    if mod(i, 10) == 0 
        % Definindo os limites para o gráfico
        k1_range_opt = unique(positions_opt(:, 1)); 
        k2_range_opt = unique(positions_opt(:, 2));
        [K1_opt, K2_opt] = meshgrid(k1_range_opt, k2_range_opt); 
        Z_opt = griddata(positions_opt(:, 1), positions_opt(:, 2), positions_opt(:, 3), K1_opt, K2_opt); 

        % Gráfico de Superfície (surf)
        figure(89);
        surf(K1_opt, K2_opt, Z_opt); % Plota a superfície
        title('Superfície das Posições das Partículas Após Otimização','FontSize', 14);
        xlabel('k1','FontSize', 12); ylabel('k2','FontSize', 12); zlabel('k3','FontSize', 12);
        shading interp; 
        colorbar; 
    end
end

%% Exibição dos resultados após a otimização
disp('Comparativo dos resultados após a otimização:');
disp(['Melhor score após: ', num2str(global_best_score_opt)]);
disp(['Melhor posição após: ', num2str(global_best_position_opt)]);
disp(['Aceleração RMS após: ', num2str(a_rms_after)]);


%% Inicialização dos parâmetros
max_iter = n_iterations;
global_best_score_opt = inf; 
global_best_position = zeros(1, 5); 
convergence = zeros(max_iter, 1); 

% Inicializando as partículas e suas posições
positions = rand(n_particles, 5); 
velocities = zeros(n_particles, 5); 

% Função de fitness 
fitness_func = @(x) sum((x - [500000, 500000, 500000, 25000, 25000]).^2);

%% Iniciando o loop de iteração
for iter = 1:max_iter
    for i = 1:n_particles
        % Atualiza a posição e a velocidade das partículas 
        velocities(i, :) = 0.5 * velocities(i, :) + 1.5 * rand(1, 5) .* (positions(i, :) - global_best_position);
        positions(i, :) = positions(i, :) + velocities(i, :);
        
        % Limita as posições das partículas para dentro de limites específicos
        positions(i, :) = max(min(positions(i, :), [1000000, 1000000, 1000000, 50000, 50000]), [0, 0, 0, 0, 0]);
        
        % Calcula o fitness da partícula
        fitness_value = fitness_func(positions(i, :));
        
        % Atualiza o melhor score global e a melhor posição global se necessário
        if fitness_value < global_best_score_opt
            global_best_score_opt = fitness_value;
            global_best_position = positions(i, :); % Armazena a posição que obteve o melhor fitness
        end
    end
    
    % Armazena o melhor fitness global a cada iteração
    convergence(iter) = global_best_score_opt;
end

%% Plot da Convergência
figure;
plot(1:max_iter, convergence, 'LineWidth', 2);
title('Convergência do PSO');
xlabel('Iterações');
ylabel('Melhor Fitness Global');
grid on;

%% Visualização das Posições das Partículas
figure;
plot3(positions(:, 1), positions(:, 2), positions(:, 3), 'o');
title('Posições das Partículas');
xlabel('k1');
ylabel('k2');
zlabel('k3');
grid on;

%% Atualização dos gráficos de comparação
figure (15)
hold on;
plot( t, As, 'b','LineWidth', 1.5, 'DisplayName', 'Antes da Otimização');
plot( t, As_opt, 'r','LineWidth', 1.3, 'DisplayName', 'Após a Otimização');
title('Comparação da Aceleração Vertical Antes e Após Otimização', 'FontSize', 14);
xlabel('t (s)', 'FontSize', 12);
ylabel('a (m/s²)', 'FontSize', 12);
legend('show');
hold off;
grid on;
set(gca, 'FontSize', 12); % Aumenta o tamanho da fonte dos eixos

%% Matriz de estado para o sistema massa-mola-amortecedor (assento)
A_assento = [0 1; -k_assento/m_assento -b_assento/m_assento];
B_assento = [0; 1];
C_assento = [1 0]; % Saída: deslocamento do assento
D_assento = 0;
a_chassi_opt = As_opt;
% Criar sistema de espaço de estados
sys_assento = ss(A_assento, B_assento, C_assento, D_assento);

% Resolver a resposta dinâmica do assento usando a aceleração do chassi como entrada
Z_assento_opt = lsim(sys_assento, m_assento * a_chassi_opt, t);
t_span=t';

% Matrizes de Estado para Cada Segmento
A_1 = [0, 1; -K_1/M_1, -B_1/M_1];  % Matriz para a cabeça
A_2 = [0, 1; -K_2/M_2, -B_2/M_2];  % Matriz para o torso superior
A_3 = [0, 1; -K_3/M_3, -B_3/M_3];  % Matriz para o torso inferior
A_4 = [0, 1; -K_4/M_4, -B_4/M_4];  % Matriz para as coxas/pélvis

B_f_1 = [0; 1/M_1];  % Força aplicada no segmento 1 (cabeça)
B_f_2 = [0; 1/M_2];  % Força aplicada no segmento 2 (torso superior)
B_f_3 = [0; 1/M_3];  % Força aplicada no segmento 3 (torso inferior)
B_f_4 = [0; 1/M_4];  % Força aplicada no segmento 4 (coxa/pélvis)

% Definindo a matriz 
C = [1, 0];  
D = 0;       

% Sistema de Estado para Cada Segmento (utilizando o comando ss)
sys_1_opt = ss(A_1, B_f_1, C, D);  % Sistema para o segmento 1 (cabeça)
sys_2_opt = ss(A_2, B_f_2, C, D);  % Sistema para o segmento 2 (torso superior)
sys_3_opt = ss(A_3, B_f_3, C, D);  % Sistema para o segmento 3 (torso inferior)
sys_4_opt = ss(A_4, B_f_4, C, D);  % Sistema para o segmento 4 (coxa/pélvis)

x1_opt = lsim(sys_1_opt, Z_assento_opt, t_span);  % Resposta para o segmento 1 (cabeça)
x2_opt = lsim(sys_2_opt, Z_assento_opt, t_span);  % Resposta para o segmento 2 (torso superior)
x3_opt = lsim(sys_3_opt, Z_assento_opt, t_span);  % Resposta para o segmento 3 (torso inferior)
x4_opt = lsim(sys_4_opt, Z_assento_opt, t_span);  % Resposta para o segmento 4 (coxa/pélvis)

%% Extrair a deflexão e a aceleração 
deflexao1_opt = x1_opt(:, 1);  % Deflexão da cabeça
deflexao2_opt = x2_opt(:, 1);  % Deflexão do torso superior
deflexao3_opt = x3_opt(:, 1);  % Deflexão do torso inferior
deflexao4_opt = x4_opt(:, 1);  % Deflexão das coxas/pélvis

figure
plot(t, deflexao1_opt, t, deflexao2_opt, t, deflexao3_opt, t, deflexao4_opt);
title('Deflexão dos Segmentos do Corpo Otimizado');
legend('Cabeça', 'Torso Superior', 'Torso Inferior', 'Coxas/Pélvis');
grid on
xlabel('t (s)');
ylabel('D (m)');

%% Encontrar o máximo deslocamento de cada parte
max_deslocamento_cabeca_opt = max(deflexao1_opt);
max_deslocamento_torso_superior_opt = max(deflexao2_opt);
max_deslocamento_torso_inferior_opt = max(deflexao3_opt);
max_deslocamento_pernas_opt = max(deflexao4_opt);

% Exibir os máximos deslocamentos
disp(['Máximo Deslocamento da Cabeça Após Otimização: ', num2str(max_deslocamento_cabeca_opt), ' m']);
disp(['Máximo Deslocamento do Torso Superior Após Otimização: ', num2str(max_deslocamento_torso_superior_opt), ' m']);
disp(['Máximo Deslocamento do Torso Inferior Após Otimização: ', num2str(max_deslocamento_torso_inferior_opt), ' m']);
disp(['Máximo Deslocamento das Pernas Após Otimização: ', num2str(max_deslocamento_pernas_opt), ' m']);


%% Lógica condicional para determinar a categoria de conforto
% Análise de Conforto para Aceleração RMS 
        if a_rms_before < 0.315
            disp('Confortável') 
        elseif a_rms_before >= 0.315 && a_rms_before  < 0.63
            disp('Um pouco desconfortável') 
        elseif  a_rms_before  >= 0.5 && a_rms_before < 1
        disp ('Relativamente desconfortável') ;
            elseif a_rms_before >= 0.8 && a_rms_before < 1.6
        disp ('Desconfortável') ;
    elseif a_rms_before >= 1.25 &&  a_rms_before < 2.5
        disp('Muito desconfortável') ;
    else
        disp ('Extremamente desconfortável') ;
    end

 % Análise de Conforto para Aceleração RMS Após Otimização
        if a_rms_after < 0.315
            disp('Confortável') 
        elseif a_rms_after  >= 0.315 && a_rms_after   < 0.63
            disp('Um pouco desconfortável') 
        elseif a_rms_after   >= 0.5 &&  a_rms_after < 1
        disp ('Relativamente desconfortável') ;
        elseif a_rms_after  >= 0.8 && a_rms_after < 1.6
        disp ('Desconfortável') ;
    elseif a_rms_after  >= 1.25 &&  a_rms_after < 2.5
        disp('Muito desconfortável') ;
    else
        disp ('Extremamente desconfortável') ;
        end   


%% Apresentação VDV
    t1=t(:);    
delta_t = diff(t1);  

% Número de pontos de amostragem
nv = length(As);

% Calcular o VDV usando a fórmula
VDV = (sum((As(1:end-1) .* delta_t) .^ 4))^(1/4);

% Exibir o VDV calculado
disp(['Valor de Dose de Vibração (VDV): ', num2str(VDV), ' m/s^2']);

t1=t(:);

% Calcular o intervalo de tempo
delta_t = diff(t1);

% Número de pontos de amostragem
nvv = length(As_opt);

% Calcular o VDV otimizado
VDV_opt = (sum((As_opt(1:end-1) .* delta_t) .^ 4))^(1/4);

% Exibir o VDV calculado
disp(['Valor de Dose de Vibração (VDV) após otimização: ', num2str(VDV_opt), ' m/s^2']);


%% Gráfico de Convergência do PSO
figure;
plot(1:max_iter, convergence, 'LineWidth', 2);
title('Convergência do PSO');
xlabel('Iterações');
ylabel('Melhor Fitness Global');
grid on;

%% Gráficos Deflexão
figure
plot(t, deflexao1, 'LineWidth', 1.5);
hold on;
plot(t, deflexao1_opt, '--', 'LineWidth', 1.5);
title('Deslocamento da Cabeça');
xlabel('t (s)');
ylabel('D (m)');
legend('Calculado', 'Otimizado');
grid on;

% Comparação do deslocamento do torso superior
figure
plot(t, deflexao2, 'LineWidth', 1.5);
hold on;
plot(t, deflexao2_opt, '--', 'LineWidth', 1.5);
title('Deslocamento do Torso Superior');
xlabel('t (s)');
ylabel('D (m)');
legend('Calculado', 'Otimizado');
grid on;

% Comparação do deslocamento do torso inferior
figure
plot(t, deflexao3, 'LineWidth', 1.5);
hold on;
plot(t, deflexao3_opt, '--', 'LineWidth', 1.5);
title('Deslocamento do Torso Inferior');
xlabel('t (s)');
ylabel('D (m)');
legend('Calculado', 'Otimizado');
grid on;

% Comparação do deslocamento das pernas
figure
plot(t, deflexao4, 'LineWidth', 1.5);
hold on;
plot(t, deflexao4_opt, '--', 'LineWidth', 1.5);
title('Deslocamento das Pernas');
xlabel('t (s)');
ylabel('D (m)');
legend('Calculado', 'Otimizado');
grid on;

%% Respostas escritas
disp('antes');
disp(k1);
disp(k2);
disp(k3);
disp(br);
disp(bB);
disp ('otimizados');
disp(k1_opt);
disp(k2_opt);
disp(k3_opt);
disp(br_opt);
disp(bB_opt);
disp('máximo angulo pitch');
disp(max(pitch_angle));
disp('n particulas');
disp(n_particles);
disp('n iterações');
disp(n_iterations);