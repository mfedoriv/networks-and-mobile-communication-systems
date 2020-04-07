clc
clear
close all

% ������ 2: ������� 3
% ------------------------ �������� ������ --------------------------------
%crc16 k=8 1+x2+x5+x6+x8+x10+x11+x12+x13+x16
g_x = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1];%crc16 k=8 x16+x13+x12+x11+x10+x8+x6+x5+x2+1
% g_x = [1, 1, 0, 1]; % 1 + X + X^3
k = 8;
% k = 4;
r = length(g_x) - 1;
n = k + r;

% ������������� ����������� � ������������ ������ ���� ��������
[H, G] = hammgen(3); 
% ������������� �������� ���������� �������� ������ ��� ����������� ������������� �� ��������
hammParam = [7, 4];
eVectors = MostSuitableEVectors(H, hammParam);

E = 1;
% N = 25000; % ���������� �������������
SNRdB = -10:2:10;
SNRdB_theor_arr = -10:10;
messages = de2bi(0:1:2^k-1);%����� ����� �� 0 �� 2^k-1 ����� ������� � ������� ���������� �������

codewords = zeros(2^k, n); % ������ ������� ����
mess_6_4 = zeros(6, 4, 2^k); % ������ ������� ����, �������� �� 6 ������ �� 4 ����
hamm_6_7 = zeros(6, 7, 2^k); % ������ ���� ������������� �� ��������
mod_6_7 = zeros(6, 7, 2^k);  % ������ �������������� ������� ����
for i = 1 : 2^k
    codewords(i,:) = CRCcoder(messages(i,:), g_x);
    mess_6_4(:, :, i) = reshape(codewords(i, :), [4, 6])';
    for j = 1 : 6
        hamm_6_7(j,:,i) = mod(mess_6_4(j,:,i)*G, 2); % ����������� �� ��������
        mod_6_7(j,:,i) = hamm_6_7(j,:,i).*2-1; % ������������� BPSK
    end
end

% ----------------- ����������� ������ �� ��� ��� BPSK --------------------
Pe_bit_theor = qfunc(sqrt(2*(10.^(SNRdB_theor_arr./10))));

% ----------------- ����������� ������ ������������� CRC ------------------
% ������� ������� (���������������):
Pe_decode_theor_asimpt = ones(1,length(SNRdB_theor_arr)).*(1/2^r);

% ����� ���. ���������� ���� � ���-�� ������� ���� �������� �����
A_i = zeros(1, n + 1); % ���������� ���� ���� i
weights = sum(codewords(:,:), 2)'; % ������ ����� ���� ������� ����
for i = 1 : n+1
    A_i(i) = sum(weights(:) == (i-1));
end
d = min(weights(:, 2:end)); % ����������� ���������� ����
% ------------- ������ ����������� ������ ������������� CRC: --------------
Pe_decode_theor_exact = zeros (1, length(SNRdB_theor_arr));
for i = 1 : length(SNRdB_theor_arr)
    for j = d : n
        Pe_decode_theor_exact(1, i) = Pe_decode_theor_exact(1, i) + A_i(j + 1) * ...
            Pe_bit_theor(i)^j * (1 - Pe_bit_theor(i))^(n - j);
    end
end

% ������������� ���� ������� ���� ��� ������� ������������� �� ��������
all_cw_for_soft = zeros(2^4, 4);
h = zeros(2^k, 7);
for i = 1:2^4
    all_cw_for_soft(i, :) = de2bi(i-1, 4);
    temp = mod(all_cw_for_soft(i, :) * G, 2);
    h(i, :) = temp(:)*2-1;
end

% -------------------------- ������������� --------------------------------

ind = 1;
Pb = zeros(1, length(SNRdB)); % ������ �� ��� � ������
Pb_synd = zeros(1, length(SNRdB)); % ������ �� ��� ��� ������������� ��������
Pb_soft = zeros(1, length(SNRdB)); % ������ �� ��� ��� ������ ������������� ��������
Pe_decode = zeros(1, length(SNRdB)); % ������ ������������� CRC ��� ������������� ��������
Pe_decode_soft = zeros(1, length(SNRdB)); % ������ ������������� CRC ��� ������ ������������� ��������
T = zeros(1, length(SNRdB)); % ���������� �����������
tests = [100 400 1500 5000 10000 10000 10000 10000 10000 10000 10000]; 
time_all = 0;
disp('Start of modelling...');
for SNR = SNRdB
    % Nb - ���-�� ������� ������, Ne_decode - ���-�� ������ �������������
    % ��� ���������� ������������� �� ��������, Ne_decode_soft - ���-��
    % ������ ������������� ��� ������ ������������� �� ��������
    Nb = 0; Nb_synd = 0; Nb_soft = 0; Ne_decode = 0; Ne_decode_soft = 0; 
    SNRi = 10^(SNR/10);
    sigma = sqrt(E/(2*SNRi));
    % Nt - ���-�� ������, Nsent - ���-�� ������������ ���������
    Nt = 0; Nsent = 0;
    tic
    for i = 1:tests(ind)
        % ��������� �������� �����
        rnd_ind = randi(2^k, 1);
        % �RC + Hamming + BPSK
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

            % ���� S == 0, � ������ ������ != 0
            if (bi2de(s_hammsoftdecoder) == 0) && (e_soft > 0)
                Ne_decode_soft = Ne_decode_soft + 1; % ������ �������������(soft) + 1
            end
            % ���� S == 0
            if bi2de(s_hammsynddecoder) == 0
                % ������ ������ != 0
                if e > 0
                    Ne_decode = Ne_decode + 1; % ������ ������������� + 1
                end
                Nsent = Nsent + 1; % ���-�� ������������ ��������� + 1
                break; % ����� �� ������� �����
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
% ------------------------------ ������� ----------------------------------

function [a_x] = CRCcoder(m_x, g_x)
% ������� ������ CRC
% ���������:
% m_x - ���������
% g_x - ����������� ���������
% ���������:
% a_x - ������� �����
    n = length(m_x) + length(g_x)-1; % len of a_x
    a_x = [zeros(1, n - length(m_x)), m_x]; % a_x * x^r (�����)
    [~, c_x] = gfdeconv(a_x, g_x); % a_x mod g_x
    a_x = gfadd(c_x, a_x);
end


function eVectors = MostSuitableEVectors(H, hammParam)
% ������� ������� �������� ���������� ������� ������ ��� �������
% ��������
% ���������:
% H - ����������� ������� ���� ��������
% hammParam - ������ � ����������� ���� ��������, � ������ ������ ����� ���
%             [7, 4]
% ���������:
% eVectors - ������� �������� ���������� �������� ������ ��� �������
%            �������� � ���������� ����, �.�. ��� ����, ����� ��������
%            �������� ���������� ������ ��� �������� s = 100(bin) = 4(dec),
%            ����� ���������� � ����� ������� ��������� �������: 
%            eVectors(4, :)

    s = 1:bi2de(ones(1, hammParam(1)-hammParam(2)));
    eVectors = zeros(length(s), hammParam(1));
    for i = s
        tmp = [];
%         e = 1:2^hammParam(1)-1
        % �������� �� ���� ��������� �������� ������
        for e = 1:bi2de(ones(1, hammParam(1)))
            % ������� ���, ��������������� ������� s == e*H'
            if i == bi2de(mod(de2bi(e, hammParam(1))*H', 2))
                tmp = [tmp e];
            end
        end
        % �������������� min ����������� ��������� ��������� ��� ������
        % ���������� ��������
        min = sum(ones(1, hammParam(1)));
        for j = 1:length(tmp)
            % ���� ��� �������� ������� ������ < min
            if sum(de2bi(tmp(j))) < min
                % ������������� ���� ���� ��� �����������
                min = sum(de2bi(tmp(j)));
                % ���������� ���� ������ ��� �������� ���������
                % ������ ������ ��� �������� �������� i
                eVectors(i, :) = de2bi(tmp(j), hammParam(1));
            end
        end
    end
end


function [decodeCW, nErrBits, v] = HammingSyndromeDecoder(noisy_mes, H, eVectors, modCW)
% ������� ����������� ������������� ��������
% ���������:
% noisy_mes - ����������� ��������� �������� �� 6 ������ �� 7 ���,
%             �.�. ��� ������ �������� 6�7
% H - ����������� ������� ���� ��������
% eVectors - �������� ���������� ������� ������ ��� ������� ��������
% modCW - �������������� ������� �����
% ����������:
% decodeCW - �������������� ������� ����� �������� 24 ����
% nErrBits - ���������� ������� ������
% v - ���������� ������� ������ ����� ���� ��������
    
    nErrBits = sum(sum(xor(noisy_mes > 0, modCW > 0)));
    demodW = double(noisy_mes(:,:)>0);
%     decode = zeros(6, 7);
    for part = 1:6
        s = bi2de(mod(demodW(part, :)*H', 2)); % ���������� ��������
        if s > 0
            demodW(part, :) = gfadd(demodW(part, :), eVectors(s, :));
        end
%         decode(part, :) = demodW(part, 1:end);
    end
    % ���-�� ������ ����� ������������� �� ��������
    v = sum(sum(xor(demodW, modCW > 0)));
    decodeCW = reshape(demodW(:, 4:end)', 1, 24);
end


function [decodeCW, v] = HammingSoftDecoder(noisy_mes, allCW, h, modCW)
% ������� ������� ������������� �� ��������
% ���������:
% noisy_mes - ����������� ��������� �������� �� 6 ������ �� 7 ���,
%             �.�. ��� ������ �������� 6�7
% allCW - ��� ��������� ��������� �� 4 ���
% h - ��� ��������� ������������������
%     ����� ����������� �� �������� � �������������, �.�. ������ ������� 2^4�6
% modCW - �������������� ������� �����
% ����������:
% decodeCW - �������������� ������� ����� �������� 24 ����
% v - ���������� ������� ������ ����� ���� ��������

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