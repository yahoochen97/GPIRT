# loading packages
import numpy as np
import pandas as pd

MODEL = ["ggum_uni_", "graded_uni_", "gpcm_uni_", "sequential_uni_",
         "DSEM_", "GPDM_", "NIRT_", "LIRT_", "DOIRT_", "GPIRT_"]

RESULT_PATH = "./log/"
MAX_SEED = 25

results = np.zeros((len(MODEL), MAX_SEED))

for SEED in range(1, MAX_SEED + 1):
    for i in range(len(MODEL)):
        loading_file = RESULT_PATH + "{}20_{}.log".format(MODEL[i],SEED)

        with open(loading_file) as f:
            data = f.readlines()
            data = data[-5]
            data = int(data.split(":")[1].replace(" ", "")[:-5])
            results[i, SEED - 1] = data

results_mu = np.round(np.mean(results, axis=1), decimals=3)
results_std = np.round(np.std(results, axis=1) / np.sqrt(MAX_SEED), decimals=3)

results = pd.DataFrame({"model": MODEL, "avg runtime": results_mu, "std runtime": results_std})
print(results)

        

