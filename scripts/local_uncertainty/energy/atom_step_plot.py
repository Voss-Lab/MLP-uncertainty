"""
This script plots the first occurrence of spikes for atoms in three bonding regimes: gas phase H (H2), Pt, and surface H (H*).
"""
import matplotlib.pyplot as plt

# Data for each bin from spike.py
# bin1, bin2, and bin3 correspond to first occurence of spikes for atoms of gas phase H, Pt, and surface H respectively
bin1 = {1188: 0, 1189: 0, 1278: 0, 1279: 0, 1187: 47, 1177: 125, 1206: 399, 1207: 555, 1159: 706, 1158: 754, 1181: 833, 1306: 836}
bin2 = {72: 451, 39: 572, 17: 617, 50: 623, 51: 626, 62: 630, 6: 778, 829: 801, 817: 802, 841: 833, 840: 835, 818: 848, 828: 867, 806: 919, 845: 925, 795: 948, 136: 958, 833: 958, 834: 964, 844: 977}
bin3 = {1047: 595, 962: 789, 1059: 832, 989: 948, 939: 953}

bins = {'$H_2$': bin1, 'Pt': bin2, 'H*': bin3}
colors = ['blue', 'green', 'red']
fig, ax = plt.subplots(figsize=(10, 8))

for idx, (bin_name, bin_data) in enumerate(bins.items()):
    ax.scatter(bin_data.values(), bin_data.keys(), color=colors[idx], label=bin_name, alpha=0.7)
ax.set_xlabel("MD steps", fontsize=25)
ax.set_ylabel("Atom ID", fontsize=25)
ax.legend(fontsize=25)
ax.set_xlim([0,1000])
ax.set_ylim([0,1400])
ax.set_yticks(range(0, max(bin1.keys()) + 100, 200))
ax.tick_params(axis='both', which='major', labelsize=20)
plt.savefig("spike.png", format='png', dpi=300)
plt.show()