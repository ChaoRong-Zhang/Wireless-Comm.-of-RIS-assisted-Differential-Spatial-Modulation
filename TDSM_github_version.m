clc;
clear;
tic;
format shortE;

Nt = 3; 
Nr = [1]; % Can be changed to an array to test different numbers of antennas, e.g., [1, 2, 4]
M = 2;    % PSK modulation order (M=2 means BPSK)
b = log2(M);
SNR = -20:2:0;
N = 64;   % Number of RIS elements

% Pre-calculate lookup tables and permutations (significantly speeds up simulation)
len1 = floor(log2(factorial(Nt))); % Number of bits for the spatial (D) part
len2 = Nt * b;                     % Number of bits for the PSK part
len = len1 + len2;                 % Total number of bits transmitted per block
blocks_per_frame = 100;            % Number of blocks transmitted per frame

% Mapping table for the spatial (D) part
tot_tables = 2^(len1); 
temp_arr = 1:Nt; 
lookup_all = flip(perms(temp_arr), 1); 
lookup = lookup_all(1:tot_tables, :); % Extract valid spatial mapping combinations

% Pre-generate candidate matrix combinations for the receiver
num_perm = len1; 
num_bits = Nt * b; 
permutations = lookup; % Receiver only compares valid spatial combinations
bin_com_dec = 0:(2^num_bits - 1); 
bin_com = de2bi(bin_com_dec, num_bits, 'left-msb'); % Generate all possible PSK bit sequences

% Matrix to record BER
berDSM = zeros(length(Nr), length(SNR));

for i_nr = 1:length(Nr)
    Nr_val = Nr(i_nr);
    
    for n = 1:length(SNR)
        numErrs = 0;
        numBits = 0;
        Sigma = 10^(-SNR(n)/10);
        
        % Modified termination condition: stop after 500 errors or 1e6 bits to ensure statistical stability
        while numErrs < 500 && numBits < 1e6 
            % Generate input data
            dataIn = randi([0 1], 1, blocks_per_frame * len); 
            dataOut_frame = zeros(1, blocks_per_frame * len); % Pre-allocate output array instead of string concatenation
            
            %% Establish Block Fading Channel (Channel remains constant within a frame)
            % H1: Tx to RIS channel (N x Nt), Rayleigh fading
            H1 = sqrt(0.5) * (randn(N, Nt) + 1i*randn(N, Nt));
            
            % H2: RIS to Rx channel (Nr x N), Nakagami-m fading (m=1.5)
            m_naka = 1.5;
            omega = 1;
            env = sqrt(gamrnd(m_naka, omega/m_naka, Nr_val, N));
            phase = exp(1i * 2 * pi * rand(Nr_val, N));
            H2 = env .* phase; 
            
            % RIS phase alignment (randomly select one Tx antenna for alignment)
            xu = randi(Nt); 
            % Assume phase alignment compensation based on the 1st Rx antenna
            G_diag = exp(-1i * angle(H2(1, :).' .* H1(:, xu))); 
            Phi = diag(G_diag); % N x N
            
            % Calculate the cascaded equivalent channel Heff (Nr x Nt)
            Heff = H2 * Phi * H1; 
            
            %% Transmit Reference Block (Differential modulation requires an initial reference block)
            St_prev = eye(Nt); % Initial reference state
            AWGN_prev = sqrt(Sigma/2) * (randn(Nr_val, Nt) + 1i*randn(Nr_val, Nt));
            Yt_prev = Heff * St_prev + AWGN_prev; % Reference signal received at the Rx (Nr x Nt)
            
            %% Block-by-Block Transmission within the Frame
            for t = 1:blocks_per_frame
                x = dataIn(len*(t-1)+1 : len*t); % Bits transmitted in the current block
                bits1 = x(1:len1);               % Bits for the spatial part
                bits2 = x(len1+1 : len1+len2);   % Bits for the PSK part
                
                % Lookup table and constellation mapping
                ind = lookup(bi2de(bits1, 'left-msb') + 1, :); 
                matrix = full(ind2vec(ind, Nt)); 
                
                Xt = matrix;
                for col = 1:Nt
                    % Extract and modulate PSK bits corresponding to the active antenna
                    psk_bits = bits2((col-1)*b + 1 : col*b);
                    symbol = pskmod(bi2de(psk_bits, 'left-msb'), 2^b);
                    Xt(:, col) = Xt(:, col) * symbol; 
                end
                
                % Differential encoding
                St = St_prev * Xt; 
                
                % Pass through the channel and add independent noise matrix
                AWGN = sqrt(Sigma/2) * (randn(Nr_val, Nt) + 1i*randn(Nr_val, Nt));
                Yt = Heff * St + AWGN; % Yt dimension is Nr x Nt
                
                %% Receiver Detection (ML Detector)
                y = Yt' * Yt_prev; % Core matrix for differential detection, dimension is Nt x Nt
                
                max_val = -Inf; % Fix the bug where the built-in 'max' function was overwritten
                final_idx = 1;
                X_t_est = zeros(Nt);
                
                % Iterate over all possible valid mapping tables and PSK combinations
                for idx = 1:(2^num_perm)
                    for j = 1:size(bin_com, 1)
                        mat1 = zeros(Nt);
                        for k = 1:Nt
                            % Get the PSK symbol for this combination
                            val = pskmod(bi2de(bin_com(j, (k-1)*b+1 : k*b), 'left-msb'), 2^b);
                            mat1(permutations(idx, k), k) = val; 
                        end
                        
                        temp_mat = y * mat1;
                        temp_val = trace(real(temp_mat)); % Differential ML decision metric
                        
                        if temp_val > max_val 
                            max_val = temp_val;
                            X_t_est = mat1;
                            final_idx = idx;
                        end
                    end
                end
                
                % Bit de-mapping
                bits_D_est = de2bi(final_idx - 1, num_perm, 'left-msb');
                bits_PSK_est = zeros(1, num_bits);
                for idx_nt = 1:Nt
                    est_sym = X_t_est(permutations(final_idx, idx_nt), idx_nt); 
                    bits_PSK_est((idx_nt-1)*b + 1 : idx_nt*b) = de2bi(pskdemod(est_sym, 2^b), b, 'left-msb');
                end
                
                % Store in the output frame
                dataOut_frame((t-1)*len + 1 : t*len) = [bits_D_est, bits_PSK_est];
                
                % Update states
                St_prev = St; 
                Yt_prev = Yt; 
            end
            
            %% Bit Error Rate (BER) Calculation
            % Compare errors for the entire frame (including the first data block)
            nErrors = biterr(dataIn, dataOut_frame);
            numErrs = numErrs + nErrors;
            numBits = numBits + length(dataIn);
        end
        
        berDSM(i_nr, n) = numErrs / numBits;
        disp(berDSM);
    end
end
toc;

%% Plotting Section
figure;
semilogy(SNR, berDSM(1, :), '-rp', 'LineWidth', 1.5, 'MarkerSize', 8); 
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('RIS-assisted Differential Spatial Modulation');
legend(sprintf('DSM N_t=%d, N_r=%d, b=%d', Nt, Nr(1), b));