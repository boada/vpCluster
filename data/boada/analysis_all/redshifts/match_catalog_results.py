import numpy as np
import pandas as pd
import sys

def main(cluster):
    results = ['results_may12', 'results_august12', 'results_may13']
    cat = pd.read_csv('./../catalogs/'+cluster+'_complete.csv')

    for result in results:
        print result
        redshifts_part = np.genfromtxt(result, names=True,
                dtype=None)

        # filter out for the cluster that we want
        mask = redshifts_part['cluster'] == cluster
        redshifts_part = redshifts_part[mask]
        redshifts_part = pd.DataFrame(redshifts_part)

        try:
            redshifts = redshifts.append(redshifts_part)
        except NameError:
            redshifts = redshifts_part

    merged = pd.merge(cat, redshifts, left_on=['tile', 'dither', 'fiber'],
            right_on=['tile','dither','fiber'], how='outer')

    try:
        del merged['index']
    except KeyError:
        #delete the first column
        del merged['Unnamed: 0']

    merged.to_csv(cluster+'_redshifts.csv', index=False)
if __name__ == "__main__":
    cluster = sys.argv[1]
    main(cluster)
