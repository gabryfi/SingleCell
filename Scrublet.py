# installazione dipendenze e pacchetti
pip install scrublet matplotlib numpy scipy scikit-learn

# importa pacchetti
import scrublet as scr
import matplotlib.pyplot as plt
import numpy as np


# Given a raw (unnormalized) UMI counts matrix counts_matrix with cells as rows and genes as columns, calculate a doublet score for each cell

import scrublet as scr
import matplotlib.pyplot as plt

scrub = scr.Scrublet(counts_matrix)

doublet_scores, predicted_doublets = scrub.scrub_doublets()

scrub.plot_histogram()
plt.show()

# score basso → probabilmente singlet
# score alto → probabilmente doublet


# rimuovi i doublets
filtered_counts_matrix = counts_matrix[~predicted_doublets]
