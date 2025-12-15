function cost = calculate_DE_cost(w_real_imag, M, Aq, PdM, alpha, weights)

% w_real_imag: Vector chứa phần thực và phần ảo của các trọng số w
% M: Số phần tử ăng-ten
% Aq: Ma trận đáp ứng mảng
% PdM: Búp sóng mong muốn (có thể là giá trị phức)
% alpha: Các chỉ số hướng cần tối ưu
% weights: Trọng số cho các điểm sai số

% 1. Tái tạo lại vector trọng số phức w từ vector thực và chuyển vị thành vector cột
w_complex = (w_real_imag(1:M) + 1j * w_real_imag(M+1:end)).';

% 2. Chuẩn hóa vector w để đảm bảo ràng buộc công suất (tương đương cₛ)
w_complex = w_complex / norm(w_complex);

% 3. Tính búp sóng phức được tạo ra tại các hướng alpha (tương đương A*w)
pattern_generated_complex = w_complex' * Aq(:, alpha);

% 4. Lấy búp sóng mong muốn tại các hướng alpha (tương đương v)
pattern_desired = PdM(alpha);

% 5. Tính sai số phức có trọng số (tương đương D*(A*w - v))
error_vector = pattern_generated_complex - pattern_desired;

% 6. Tính chuẩn bậc 2 bình phương của sai số có trọng số ||D(Aw - v)||₂²
cost = sum(weights .* abs(error_vector).^2);

end
