# Plotting empirical probabilities

Permutation and spin annotations are recalculated from their null draws and observed statistics using `mean(abs(null) >= abs(observed))`, with inclusive ties and no +1 correction. A stored p field is not used for the distribution plots. Joint spin/permutation panels retain the maximum of their separately recalculated probabilities. Bootstrap panels retain their specified bootstrap definition.

Replacement random-order tests use the lower tail, `mean(null <= observed)`, without a +1 correction, followed by Holm adjustment across the 18 comparisons.
