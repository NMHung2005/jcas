function [best_w_vec] = optimize_DE(objective_func, M, de_params)
% HÀM THỰC THI THUẬT TOÁN DIFFERENTIAL EVOLUTION

% --- Trích xuất tham số DE ---
pop_size = de_params.pop_size; % Kích thước quần thể
F = de_params.F;               % Hệ số đột biến
CR = de_params.CR;             % Tỷ lệ lai ghép
max_gens = de_params.max_gens; % Số thế hệ tối đa
dim = 2 * M;                   % Số chiều (phần thực + phần ảo)

% --- Định nghĩa không gian tìm kiếm ---
lower_bound = -1 * ones(1, dim);
upper_bound =  1 * ones(1, dim);

% 1. KHỞI TẠO QUẦN THỂ
% Tạo ngẫu nhiên các vector trọng số trong không gian tìm kiếm
population = lower_bound + rand(pop_size, dim) .* (upper_bound - lower_bound);

% Đánh giá độ lỗi ban đầu cho từng cá thể
costs = zeros(pop_size, 1);
for i = 1:pop_size
    costs(i) = objective_func(population(i, :));
end

% 2. VÒNG LẶP TỐI ƯU
fprintf('Bat dau toi uu bang DE...\n');
for gen = 1:max_gens
    for i = 1:pop_size
        % --- Đột biến (Mutation) ---
        % Chọn 3 cá thể ngẫu nhiên r1, r2, r3 (khác i)
        idxs = randperm(pop_size);
        idxs(idxs == i) = [];
        r1 = idxs(1); r2 = idxs(2); r3 = idxs(3);
        
        mutant_vector = population(r1, :) + F * (population(r2, :) - population(r3, :));
        
        % Xử lý nếu vượt ngoài không gian tìm kiếm
        mutant_vector = max(mutant_vector, lower_bound);
        mutant_vector = min(mutant_vector, upper_bound);
        
        % --- Lai ghép (Crossover) ---
        trial_vector = population(i, :);
        j_rand = randi(dim); % Đảm bảo ít nhất 1 chiều thay đổi
        for j = 1:dim
            if rand() < CR || j == j_rand
                trial_vector(j) = mutant_vector(j);
            else
                trial_vector(j) = population(i, j);
            end
        end
        
        % --- Lựa chọn (Selection) ---
        trial_cost = objective_func(trial_vector);
        
        % Nếu cá thể mới tốt hơn cá thể cũ, thay thế nó
        if trial_cost < costs(i)
            population(i, :) = trial_vector;
            costs(i) = trial_cost;
        end
    end
    
    % Hiển thị tiến trình (tùy chọn)
    if mod(gen, 5) == 0
        fprintf('The he %d, Sai so tot nhat: %f\n', gen, min(costs));
    end
end

% 3. TRẢ VỀ KẾT QUẢ TỐT NHẤT
[~, best_idx] = min(costs);
best_w_vec = population(best_idx, :);
fprintf('DE hoan thanh!\n');

end
