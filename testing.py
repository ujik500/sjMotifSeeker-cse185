import logomaker
import pandas as pd
import matplotlib.pyplot as plt

trimmed_pfm = [[3, 6, 0, 0], [9, 0, 1, 0], [0, 10, 0, 0], [0, 10, 0, 0], [10, 0, 0, 0], [0, 10, 0, 0], [8, 0, 0, 0], [0, 6, 0, 0]]
for row in trimmed_pfm:
    total = sum(row)
    for i in range(len(row)):
        row[i] = row[i] / total
print(trimmed_pfm)
''' 
trimmed_pfm = pd.DataFrame(trimmed_pfm, dtype = int)
trimmed_pfm.columns = ["A", "C", "G", "T"]
print(trimmed_pfm)
ss_logo = logomaker.Logo(trimmed_pfm,
                         width=.4,
                         vpad=.05,
                         fade_probabilities=False,
                         stack_order='small_on_top',
                         color_scheme='classic')

# style using Logo methods
ss_logo.style_spines(spines=['left', 'right'], visible=False)

# style using Axes methods
ss_logo.ax.set_xticks(range(len(trimmed_pfm)))
ss_logo.ax.set_xticklabels('%+d'%x for x in [-3, -2, -1, 1, 2, 3, 4, 5])
ss_logo.ax.set_yticks([0, .5, 1])
ss_logo.ax.axvline(2.5, color='k', linewidth=1, linestyle=':')
ss_logo.ax.set_ylabel('probability')

plt.savefig('logo.pdf')
plt.show()
'''