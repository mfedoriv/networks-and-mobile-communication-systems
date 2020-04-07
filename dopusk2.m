clc
clear
close all

% Допуск 2: вариант 3
% ------------------------ Исходные данные --------------------------------
%crc16 k=8 1+x2+x5+x6+x8+x10+x11+x12+x13+x16
g_x = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1];%crc16 k=8 x16+x13+x12+x11+x10+x8+x6+x5+x2+1
% g_x = [1, 1, 0, 1]; % 1 + X + X^3
k = 8;
% k = 4;
r = length(g_x) - 1;
n = k + r;

% генерирование проверочной и попрождающей матриц кода Хэмминга
[H, G] = hammgen(3); 
% генерирование наиболее подходящих векторов ошибки для синдромного декодирования по Хэммингу
hammParam = [7, 4];
eVectors = MostSuitableEVectors(H, hammParam);

E = 1;
% N = 25000; % количество экспериментов
SNRdB = -10:2:10;
SNRdB_theor_arr = -10:10;
messages = de2bi(0:1:2^k-1);%генер сообщ от 0 до 2^k-1 слева направо в порядке увелечения степени

codewords = zeros(2^k, n); % массив кодовых слов
mess_6_4 = zeros(6, 4, 2^k); % массив кодовых слов, разбитых на 6 частей по 4 бита
hamm_6_7 = zeros(6, 7, 2^k); % массив слов закодированых по Хэммингу
mod_6_7 = zeros(6, 7, 2^k);  % массив модулированных кодовых слов
for i = 1 : 2^k
    codewords(i,:) = CRCcoder(messages(i,:), g_x);
    mess_6_4(:, :, i) = reshape(codewords(i, :), [4, 6])';
    for j = 1 : 6
        hamm_6_7(j,:,i) = mod(mess_6_4(j,:,i)*G, 2); % кодирование по Хэммингу
        mod_6_7(j,:,i) = hamm_6_7(j,:,i).*2-1; % модулирование BPSK
    end
end

% ----------------- Вероятность ошибки на бит для BPSK --------------------
Pe_bit_theor = qfunc(sqrt(2*(10.^(SNRdB_theor_arr./10))));

% ----------------- Вероятность ошибки декодирования CRC ------------------
% Верхняя граница (асимптотическая):
Pe_decode_theor_asimpt = ones(1,length(SNRdB_theor_arr)).*(1/2^r);

% Поиск мин. расстояния кода и кол-ва кодовых слов заданной длины
A_i = zeros(1, n + 1); % количество слов веса i
weights = sum(codewords(:,:), 2)'; % массив весов всех кодовых слов
for i = 1 : n+1
    A_i(i) = sum(weights(:) == (i-1));
end
d = min(weights(:, 2:end)); % минимальное расстояние кода
% ------------- Точная вероятность ошибки декодирования CRC: --------------
Pe_decode_theor_exact = zeros (1, length(SNRdB_theor_arr));
for i = 1 : length(SNRdB_theor_arr)
    for j = d : n
        Pe_decode_theor_exact(1, i) = Pe_decode_theor_exact(1, i) + A_i(j + 1) * ...
            Pe_bit_theor(i)^j * (1 - Pe_bit_theor(i))^(n - j);
    end
end

% Генерирование всех кодовых слов для мягкого декодирования по Хэммингу
all_cw_for_soft = zeros(2^4, 4);
h = zeros(2^k, 7);
for i = 1:2^4
    all_cw_for_soft(i, :) = de2bi(i-1, 4);
    temp = mod(all_cw_for_soft(i, :) * G, 2);
    h(i, :) = temp(:)*2-1;
end

% -------------------------- Моделирование --------------------------------

ind = 1;
Pb = zeros(1, length(SNRdB)); % ошибка на бит в канале
Pb_synd = zeros(1, length(SNRdB)); % ошибка на бит при декодировании Хэмминга
Pb_soft = zeros(1, length(SNRdB)); % ошибка на бит при мягком декодировании Хэмминга
Pe_decode = zeros(1, length(SNRdB)); % ошибка декодирования CRC при декодировании Хэмминга
Pe_decode_soft = zeros(1, length(SNRdB)); % ошибка декодирования CRC при мягком декодировании Хэмминга
T = zeros(1, length(SNRdB)); % пропускная способность
tests = [100 400 1500 5000 10000 10000 10000 10000 10000 10000 10000]; 
time_all = 0;
disp('Start of modelling...');
for SNR = SNRdB
    % Nb - кол-во битовых ошибок, Ne_decode - кол-во ошибок декодирования
    % при синдромном декодировании по Хэммингу, Ne_decode_soft - кол-во
    % ошибок декодирования при мягком декодировании по Хэммингу
    Nb = 0; Nb_synd = 0; Nb_soft = 0; Ne_decode = 0; Ne_decode_soft = 0; 
    SNRi = 10^(SNR/10);
    sigma = sqrt(E/(2*SNRi));
    % Nt - кол-во тестов, Nsent - кол-во отправленных сообщений
    Nt = 0; Nsent = 0;
    tic
    for i = 1:tests(ind)
        % генерация кодового слова
        rnd_ind = randi(2^k, 1);
        % СRC + Hamming + BPSK
        mod_mes = mod_6_7(:, :, rnd_ind);
        % CRC
        a = codewords(rnd_ind, :);
        while 1
            Nt = Nt + 1;
            % AVGN
            noisy_mes = mod_mes + sigma*randn(6, 7);
            % BPSK decode
            demod_mes = noisy_mes > 0;
            % Hamming decode
            [hamm_synd_decode_mes, nErrBits, v_synd] = HammingSyndromeDecoder(noisy_mes, H, eVectors, mod_mes);
            Nb_synd = Nb_synd + v_synd;
            Nb = Nb + nErrBits;
            % Hamming soft decode
            [hamm_soft_decode_mes, v_soft] = HammingSoftDecoder(noisy_mes, all_cw_for_soft, h, mod_mes);
            Nb_soft = Nb_soft + v_soft;
            % CRC decode
            [~, s_hammsynddecoder] = gfdeconv(hamm_synd_decode_mes, g_x);
            e = sum(xor(a, hamm_synd_decode_mes)); % e-vect with Hamm decode
            [~, s_hammsoftdecoder] = gfdeconv(hamm_soft_decode_mes, g_x);
            e_soft = sum(xor(a, hamm_soft_decode_mes)); % e-vect with SoftHamm decode

            % если S == 0, а вектор ошибок != 0
            if (bi2de(s_hammsoftdecoder) == 0) && (e_soft > 0)
                Ne_decode_soft = Ne_decode_soft + 1; % ошибка декодирования(soft) + 1
            end
            % если S == 0
            if bi2de(s_hammsynddecoder) == 0
                % вектор ошибок != 0
                if e > 0
                    Ne_decode = Ne_decode + 1; % ошибка декодирования + 1
                end
                Nsent = Nsent + 1; % кол-во отправленных сообщений + 1
                break; % выход из вечного цикла
            end
        end
        
    end
    
    Pb(ind) = Nb/(Nt * 42);
    Pb_synd(ind) = Nb_synd/(Nt * 42);
    Pb_soft(ind) = Nb_soft/(Nt * 42);
    Pe_decode(ind) = Ne_decode/Nt;
    Pe_decode_soft(ind) = Ne_decode_soft/Nt;
    T(ind) = (k * Nsent)/(42 * Nt);
    
    time = toc;
    fprintf('Time of modelling %d SNRdb: %f seconds\n', SNR, time)
    time_all = time_all + time;
    
    ind = ind + 1;
end
time_all = datestr(seconds(time_all),'HH:MM:SS');
fprintf('All time of modelling: %s', time_all);

% ------------------------------ Plotting ---------------------------------
figure(1);
semilogy(snr_theor_arr, Pe_decode_Theor, 'g', ...
         snr_theor_arr, P_ed, 'r', ...
         snr_arr, Pe_decode, 'bo', ...
         snr_arr, Pe_decode_soft, 'ko');
legend('Ped upper boound', 'Ped theor', 'Ped synd', 'Ped soft');
grid on;
xlabel('E/N_0, dB')

figure(2);
semilogy(snr_theor_arr, Pb_theor, ...
         snr_arr, Pb, 'o', ...
         snr_arr, Pb_synd, '*', ...
         snr_arr, Pb_soft, 'x')
legend('Pb_theor', 'Pb', 'Pb_synd', 'Pb_soft')
grid on;
xlabel('E/N_0, dB')

figure(3);
semilogy(snr_arr, T);
% ------------------------------ Функции ----------------------------------

function [a_x] = CRCcoder(m_x, g_x)
% Функция кодера CRC
% Аргументы:
% m_x - сообщение
% g_x - порождающий многочлен
% Результат:
% a_x - кодовое слово
    n = length(m_x) + length(g_x)-1; % len of a_x
    a_x = [zeros(1, n - length(m_x)), m_x]; % a_x * x^r (сдвиг)
    [~, c_x] = gfdeconv(a_x, g_x); % a_x mod g_x
    a_x = gfadd(c_x, a_x);
end


function eVectors = MostSuitableEVectors(H, hammParam)
% Функция находит наиболее подходящие вектора ошибок для каждого
% синдрома
% Аргументы:
% H - проверочная матрица кода Хэмминга
% hammParam - массив с параметрами кода Хэмминга, в данном случае имеет вид
%             [7, 4]
% Результат:
% eVectors - таблица наиболее подходящих векторов ошибки для каждого
%            синдрома в десятичном виде, т.е. для того, чтобы получить
%            наиболее подходящий вектор для синдрома s = 100(bin) = 4(dec),
%            нужно обратиться к этому массиву следующим образом: 
%            eVectors(4, :)

    s = 1:bi2de(ones(1, hammParam(1)-hammParam(2)));
    eVectors = zeros(length(s), hammParam(1));
    for i = s
        tmp = [];
%         e = 1:2^hammParam(1)-1
        % проходим по всем возможным векторам ошибок
        for e = 1:bi2de(ones(1, hammParam(1)))
            % находим все, удовлетворяющие условию s == e*H'
            if i == bi2de(mod(de2bi(e, hammParam(1))*H', 2))
                tmp = [tmp e];
            end
        end
        % инициализируем min максимально возможным значением для поиска
        % настоящего минимума
        min = sum(ones(1, hammParam(1)));
        for j = 1:length(tmp)
            % если вес текущего вектора ошибок < min
            if sum(de2bi(tmp(j))) < min
                % устанавливаем этот весь как минимальный
                min = sum(de2bi(tmp(j)));
                % записываем этот вектор как наиболее вероятный
                % вектор ошибок для текущего синдрома i
                eVectors(i, :) = de2bi(tmp(j), hammParam(1));
            end
        end
    end
end


function [decodeCW, nErrBits, v] = HammingSyndromeDecoder(noisy_mes, H, eVectors, modCW)
% Функция синдромного декодирования Хэмминга
% Аргументы:
% noisy_mes - зашумленное сообщение разбитое на 6 частей по 7 бит,
%             т.е. это массив размером 6х7
% H - проверочная матрица кода Хэмминга
% eVectors - наиболее подходящие вектора ошибок для каждого синдрома
% modCW - модулированное кодовое слово
% Результаты:
% decodeCW - декодированное кодовое слово размером 24 бита
% nErrBits - количество битовых ошибок
% v - количество битовых ошибок после кода Хэмминга
    
    nErrBits = sum(sum(xor(noisy_mes > 0, modCW > 0)));
    demodW = double(noisy_mes(:,:)>0);
%     decode = zeros(6, 7);
    for part = 1:6
        s = bi2de(mod(demodW(part, :)*H', 2)); % вычисление синдрома
        if s > 0
            demodW(part, :) = gfadd(demodW(part, :), eVectors(s, :));
        end
%         decode(part, :) = demodW(part, 1:end);
    end
    % кол-во ошибок после декодирования по Хэммингу
    v = sum(sum(xor(demodW, modCW > 0)));
    decodeCW = reshape(demodW(:, 4:end)', 1, 24);
end


function [decodeCW, v] = HammingSoftDecoder(noisy_mes, allCW, h, modCW)
% Функция мягкого декодирования по Хэммингу
% Аргументы:
% noisy_mes - зашумленное сообщение разбитое на 6 частей по 7 бит,
%             т.е. это массив размером 6х7
% allCW - все возможные сообщения по 4 бит
% h - все возможные последовательности
%     после кодирования по Хэммингу и модулирования, т.е. массив размера 2^4х6
% modCW - модулированное кодовое слово
% Результаты:
% decodeCW - декодированное кодовое слово размером 24 бита
% v - количество битовых ошибок после кода Хэмминга

    min_d = [100 100 100 100 100 100];
    decodeCW = zeros(6, 4);
    cor = zeros(6, 7);
    for part = 1:6
        for i = 1:2^4
            min = sqrt(sum((noisy_mes(part, :) - h(i, :)).^2));
            if min < min_d(part)
                min_d(part) = min;
                decodeCW(part, :) = allCW(i, :);
                cor(part, :) = h(i, :);
            end
        end
    end
    v = sum(sum(xor(cor > 0, modCW > 0)));
    decodeCW = reshape(decodeCW', 1, 24);
end