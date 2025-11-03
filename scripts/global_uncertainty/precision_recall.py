import numpy as np
import matplotlib.pyplot as plt

# Load data
sigmaE = np.loadtxt('std')  # predicted per-atom uncertainty
error = np.loadtxt('err')   # actual per-atom error

# Define high-error structures by percentile (e.g., top 20% of errors)
error_percentile = 80
error_threshold = np.percentile(error, error_percentile)
print(f"High-error threshold (top {100-error_percentile}% of errors) = {error_threshold:.4f} eV/atom")

# Define multiple sigmaE thresholds (percentiles of sigmaE)
sigma_percentiles = np.arange(50, 100, 5)  # 50th, 55th, ..., 95th
precision_list = []
recall_list = []
threshold_list = []

for p in sigma_percentiles:
    sigma_threshold = np.percentile(sigmaE, p)
    selected = sigmaE >= sigma_threshold       # predicted high-uncertainty
    high_error = error >= error_threshold      # true high-error structures

    TP = np.sum(selected & high_error)
    FP = np.sum(selected & ~high_error)
    FN = np.sum(~selected & high_error)

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0

    precision_list.append(precision)
    recall_list.append(recall)
    threshold_list.append(sigma_threshold)

    print(f"σE percentile {p} -> threshold {sigma_threshold:.4f} eV/atom: "
          f"TP={TP}, FP={FP}, FN={FN}, Precision={precision:.3f}, Recall={recall:.3f}")

# Optional: plot Precision-Recall curve
plt.figure(figsize=(7,5))
plt.plot(recall_list, precision_list, marker='o')
for i, p in enumerate(sigma_percentiles):
    plt.text(recall_list[i]+0.005, precision_list[i], f'{p}th', fontsize=10)
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title(f'Precision-Recall Curve (High-error = top {100-error_percentile}% errors)')
plt.grid(True)
plt.ylim(0,1)
plt.xlim(0,1)
plt.show()

