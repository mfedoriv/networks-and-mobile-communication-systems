clc
clear
close all

% Допуск 2: вариант 3
% ------------------------ Исходные данные --------------------------------
%crc16 k=8 1+x2+x5+x6+x8+x10+x11+x12+x13+x16
% g_x = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1];%crc16 k=8 x16+x13+x12+x11+x10+x8+x6+x5+x2+1
g_x = [1, 1, 0, 1]; % 1 + X + X^3
% k = 8;
k = 4;
r = length(g_x) - 1;
n = k + r;
% index = 1;
E = 1;
N = 25000; % количество экспериментов
% SNRdB = -30:2:0;
SNRdB = -20:2:0;
% SNRdB_theor_arr = -30:0;
SNRdB_theor_arr = -20:0;
messages = de2bi(0:1:2^k-1);%генер сообщ от 0 до 2^k-1 слева направо в порядке увелечения степени

% Генерация кодовых слов и модуляция BPSK
codewords = zeros(2^k, n);
% codewords_bpsk = zeros(2^k, n);
for i = 1 : 2^k
    codewords(i,:) = CRCcoder(messages(i,:), g_x);
%     codewords_bpsk(i,:) = codewords(i,:).*(-2) + 1;
end
% Поиск мин. расстояния кода и кол-ва кодовых слов заданной длины
A_i = zeros(1, n + 1); % количество слов веса i
weights = sum(codewords(:,:), 2)'; % массив весов всех кодовых слов
for i = 1 : n+1
    A_i(i) = sum(weights(:) == (i-1));
end
d = min(weights(:, 2:end)); % минимальное расстояние кода
% ----------------- Вероятность ошибки на бит для BPSK --------------------
Pe_bit_theor = qfunc(sqrt(2*(10.^(SNRdB_theor_arr./10))));
% ----------------- Вероятность ошибки декодирования CRC ------------------
% Верхняя граница (асимптотическая):
Pe_decode_theor_asimpt = ones(1,length(SNRdB_theor_arr)).*(1/2^r);
% Точная вероятность ошибки декодирования CRC:
Pe_decode_theor_exact = zeros (1, length(SNRdB_theor_arr));
for i = 1 : length(SNRdB_theor_arr)
    for j = d : n
        Pe_decode_theor_exact(1, i) = Pe_decode_theor_exact(1, i) + A_i(j + 1) * ...
            Pe_bit_theor(i)^j * (1 - Pe_bit_theor(i))^(n - j);
    end
end

% ----------------------- Моделирование ----------------------------------
Pe_bit = zeros (1, length(SNRdB));
Pe_decode = zeros (1, length(SNRdB));
time_all = 0;
disp('Start of modelling...');
for i = 1 : length(SNRdB)
    tic
    
    SNR = 10.^(SNRdB(i)/10);
    sigma = sqrt(E / (2*SNR));
    
    [Pe_bit(1,i), Pe_decode(1,i)] = modelling(k, g_x, codewords, sigma, N);
    
    time = toc;
    fprintf('Time of modelling %d SNRdb: %f seconds\n', SNRdB(i), time)
    time_all = time_all + time;
end
time_all = datestr(seconds(time_all),'HH:MM:SS');
fprintf('All time of modelling: %s', time_all);


% ------------------------------ Plotting ---------------------------------
f_Ped = figure(1);
semilogy(SNRdB, Pe_decode, 'b.-', ...
    SNRdB_theor_arr, Pe_decode_theor_asimpt, 'r-', ...
    SNRdB_theor_arr, Pe_decode_theor_exact, 'k-');
legend('Ped', 'Ped theor asimpt', 'Ped theor exact');
grid on;
xlabel('E/N_0, dB')

f_Pe_bit = figure(2);
semilogy(SNRdB, Pe_bit, 'b.-', ...
    SNRdB_theor_arr, Pe_bit_theor, 'r-');
legend('Pe bit', 'Pe bit theor');
grid on;
xlabel('E/N_0, dB')
% ------------------------------ Функции ----------------------------------
function [a_x] = CRCcoder(m_x, g_x)
    n = length(m_x) + length(g_x)-1; % len of a_x
    a_x = [zeros(1, n - length(m_x)), m_x]; % a_x * x^r (сдвиг)
    [~, c_x] = gfdeconv(a_x, g_x); % a_x mod g_x
    a_x = gfadd(c_x, a_x);
end

function [Pe_bit, Pe_decode] = modelling(k, g_x, codewords, sigma, N)
    r = length(g_x) - 1;
    n = k + r;
    num_of_mess = 2^k;

    Ne_bit = 0;
    Ne_decode = 0;
    
    for i = 1 : N
        % Источник
        index = randi(num_of_mess);
        % CRC coding
        a_x = codewords(index, :);
        % BPSK modulation
        a_x_bpsk = a_x.*(-2) + 1;
        % AWGN
        b =  a_x_bpsk + sigma * randn(1, n);
        
        % BPSK demodulation
        a_x_ = b < 0;
        % CRC decoding
        e = xor(a_x, a_x_);
        Ne_bit = Ne_bit + nnz(e);
        [~, s] = gfdeconv(double(a_x_), g_x);
        if (bi2de(s) == 0) && (nnz(e) ~= 0)
            Ne_decode = Ne_decode + 1;
        end
    end
    Pe_bit = Ne_bit / (N * n);
    Pe_decode = Ne_decode / N;
end