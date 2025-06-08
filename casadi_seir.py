import casadi as ca
import pandas as pd

# --- Đọc dữ liệu từ tệp Excel ---
file_path = 'C:\\Users\\Administrator\\Downloads\\data-1.xlsx'
sheet_name = 'covid19_provinces_vn'
df = pd.read_excel(file_path, sheet_name=sheet_name)

# Chọn các cột cần thiết
selected_columns = [
    "Province", "N", "S", "I", "E", "R", "Gamma", "Alpha", "Beta", "Nu", "Mu",
    "Provincial_Hospital",  "P+ICP", "Maternity_home", "C.H.C", "Vaccine_min"
]
data_filtered = df[selected_columns]

# Tạo dictionary từ dữ liệu
cities = data_filtered["Province"].tolist()
N0 = dict(zip(cities, data_filtered["N"]))
S0 = dict(zip(cities, data_filtered["S"]))
I0 = dict(zip(cities, data_filtered["I"]))
E0 = dict(zip(cities, data_filtered["E"]))
R0 = dict(zip(cities, data_filtered["R"]))
gamma_values = dict(zip(cities, data_filtered["Gamma"]))
alpha_values = dict(zip(cities, data_filtered["Alpha"]))
beta_values = dict(zip(cities, data_filtered["Beta"]))
nu_values = dict(zip(cities, data_filtered["Nu"]))
mu_values = dict(zip(cities, data_filtered["Mu"]))
Vaccine_min = dict(zip(cities, data_filtered["Vaccine_min"]))

# Tính toán năng lực phân phối vaccine C[c]
C = {}
for index, row in data_filtered.iterrows():
    C[row["Province"]] = (
        5000 * row['Provincial_Hospital'] +
        1500 * row['P+ICP'] +
        1000 * row['Maternity_home'] +
        500 * row['C.H.C']
    )

# --- Tham số ---
T = 30  # Số ngày
omega = 0.75  # Hiệu suất vaccine
V_total = 20000000  # Tổng số vaccine có thể phân phối

# --- Biến quyết định ---
opti = ca.Opti()
S = {c: opti.variable(T + 1) for c in cities}
E = {c: opti.variable(T + 1) for c in cities}
I = {c: opti.variable(T + 1) for c in cities}
R = {c: opti.variable(T + 1) for c in cities}
v = {c: opti.variable(T + 1) for c in cities}

# Tạo một dictionary để lưu dân số tại mỗi thời điểm t
N = {c: [0] * (T + 1) for c in cities}
# Gán giá trị dân số ban đầu
for c in cities:
    N[c][0] = N0[c]

# --- Hàm mục tiêu ---
objective = 0
for c in cities:
    objective += I[c][T]  # Tối thiểu hóa số ca nhiễm tại thời điểm cuối cùng
opti.minimize(objective)

# --- Ràng buộc ---
for c in cities:
    # Điều kiện ban đầu
    opti.subject_to(S[c][0] == S0[c])
    opti.subject_to(E[c][0] == E0[c])
    opti.subject_to(I[c][0] == I0[c])
    opti.subject_to(R[c][0] == R0[c])
    opti.subject_to(v[c][0] == 0)
    
    # Phương trình SEIR
    for t in range(T):
        N[c][t + 1] = N[c][t] * (1 + nu_values[c] - mu_values[c])
        opti.subject_to(S[c][t+1] == nu_values[c] * N[c][t] + S[c][t] - beta_values[c] * S[c][t] * I[c][t] / N[c][t] - omega * v[c][t+1] - mu_values[c] * S[c][t])
        opti.subject_to(E[c][t+1] == E[c][t] + beta_values[c] * S[c][t] * I[c][t] / N[c][t] - (alpha_values[c] + mu_values[c]) * E[c][t])
        opti.subject_to(I[c][t+1] == I[c][t] + alpha_values[c] * E[c][t] - (gamma_values[c] + mu_values[c]) * I[c][t])
        opti.subject_to(R[c][t+1] == R[c][t] + gamma_values[c] * I[c][t] + omega * v[c][t+1] - mu_values[c] * R[c][t])
        opti.subject_to(S[c][t+1] + E[c][t+1] + I[c][t+1] + R[c][t+1] == N[c][t+1])  # Bảo toàn dân số
        
        opti.subject_to(S[c][t+1] >= 0)
        opti.subject_to(E[c][t+1] >= 0)
        opti.subject_to(I[c][t+1] >= 0)
        opti.subject_to(R[c][t+1] >= 0)

        # Ràng buộc vaccine
        opti.subject_to(v[c][t+1] >= 0)
        opti.subject_to(v[c][t+1] <= C[c])  # Thay đổi thành ràng buộc C[c]

# Ràng buộc tổng lượng vaccine phải đạt tối thiểu Vaccine_min
for c in cities:
    opti.subject_to(ca.sum1(v[c]) >= Vaccine_min[c])

# Tổng lượng vaccine không vượt quá ngân sách
opti.subject_to(sum(ca.sum1(v[c]) for c in cities) <= V_total)

# --- Giải quyết bài toán ---
opti.solver("ipopt")  # Sử dụng solver IPOPT
sol = opti.solve()

# --- Kết quả ---
summary_results = []
daily_vaccine_results = []

for c in cities:
    total_vaccine = sum(sol.value(v[c][t]) for t in range(T))  # Tổng số vaccine cho thành phố
    summary_results.append({
        "Province": c,
        "Infected_Final": sol.value(I[c][T]),
        "Total_Vaccine_Distributed": total_vaccine,
    })
    daily_row = {"Province": c}
    for t in range(T+1):
        daily_row[f"Day_{t}"] = sol.value(v[c][t])
    daily_vaccine_results.append(daily_row)

# Xuất kết quả ra file Excel
output_file = 'C:\\Users\\Administrator\\Downloads\\SEIR_results.xlsx'

with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
    summary_df = pd.DataFrame(summary_results)
    daily_vaccine_df = pd.DataFrame(daily_vaccine_results)
    summary_df.to_excel(writer, sheet_name='Summary', index=False)
    daily_vaccine_df.to_excel(writer, sheet_name='Daily_Vaccine', index=False)

print(f"Kết quả đã được lưu tại: {output_file}")
