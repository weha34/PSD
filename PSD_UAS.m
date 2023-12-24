load('sinyal_diskrit3')

% Inisialisasi panjang sinyal dan sampling rate (fs)
fs = 100000; % Sampling rate
N = length(xn);

% Waktu sampling
t = (0:N-1) / fs;

% Plot sinyal asli (Figure 1)
figure(1);
plot(t, xn);
title('Sinyal Asli');
xlabel('Waktu (s)');
ylabel('Amplitudo');

% Hitung FFT dari sinyal xn
X = fft(xn);

% Hitung magnitudo spektrum frekuensi
mag_X = abs(X);

% Buat sumbu frekuensi
f = (0:N-1) * fs / N;

% Plot spektrum frekuensi (Figure 2)
figure(2);
stem(f, mag_X);
title('Spektrum Frekuensi');
xlabel('Frekuensi (Hz)');
ylabel('Magnitudo');

% Identifikasi frekuensi dominan
threshold = 0.5 * max(mag_X); % Ambang batas untuk frekuensi dominan
dominant_indices = find(mag_X > threshold); % Indeks frekuensi dominan

% Plot frekuensi dominan dengan warna merah
hold on;
plot(f(dominant_indices), mag_X(dominant_indices), 'ro');
hold off;

% Print indeks, frekuensi, dan magnitudo frekuensi dominan
fprintf('Frekuensi Dominan:\n');
for i = 1:length(dominant_indices)
    fprintf('Indeks: %d, Frekuensi: %.2f Hz, Magnitudo: %.2f\n', dominant_indices(i), f(dominant_indices(i)), mag_X(dominant_indices(i)));
end

% Parameter desain filter
fs = 100000; % Frekuensi sampling
f_low = 1; % Frekuensi cut-off rendah (Hz)
f_high = 4; % Frekuensi cut-off tinggi (Hz)

% Panjang filter yang berbeda
filter_lengths = [6, 12, 22, 32];

% Rancang dan plot filter untuk setiap panjang
for i = 1:length(filter_lengths)
    % Desain filter
    filter_order = filter_lengths(i) - 1;
    fc = [f_low, f_high] / (fs / 2);
    b = fir1(filter_order, fc, 'bandpass');
    
    % Plot respon frekuensi filter
    figure;
    freqz(b, 1, 1024, fs);
    title(sprintf('Frekuensi Respon Filter (Panjang %d)', filter_lengths(i)));
end

% Misalkan sinyal diskrit yang ingin difilter adalah xn
% xn merupakan sinyal yang telah dimuat sebelumnya

% Gunakan filter yang telah dirancang sebelumnya (ganti dengan koefisien filter yang sesuai)
filtered_signals = cell(length(filter_lengths), 1);

for i = 1:length(filter_lengths)
    filtered_signals{i} = filter(b, 1, xn);
end

% Plot hasil keluaran filter dalam domain waktu
figure;
for i = 1:length(filter_lengths)
    subplot(length(filter_lengths), 1, i);
    plot(filtered_signals{i});
    title(sprintf('Keluaran Filter (Panjang %d)', filter_lengths(i)));
end

% Plot hasil keluaran filter dalam domain frekuensi
figure;
for i = 1:length(filter_lengths)
    subplot(length(filter_lengths), 1, i);
    freqz(filtered_signals{i}, 1, 1024, fs);
    title(sprintf('Frekuensi Respon Filtered Signal (Panjang %d)', filter_lengths(i)));
end

nn = 1:3;
nT = 0:1/fs:0.05;
omega_l = 0.1;
omega_h = 0.3;
hn0 = omega_h - omega_l / pi;
hnkanan = (sin(omega_h * pi * nn) ./ (nn * pi)) - (sin(omega_l * pi * nn) ./ (nn * pi));
hnkiri = fliplr(hnkanan);

koef = [hnkiri, hn0, hnkanan];
N = length(nT);
disp(koef);
