from __future__ import print_function
import pandas as pd
from glob import glob
import numpy as np
from itertools import combinations
from sklearn.cross_validation import train_test_split
from sklearn.grid_search import GridSearchCV, RandomizedSearchCV
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from scipy.stats import randint as sp_randint
import h5py as hdf

# get the files
files = glob('../redshifts/*_redshifts.csv')

# load the frames into an array
frames = [pd.read_csv(f) for f in files]

# concat all of the frames together
result = pd.concat(frames, ignore_index=True)

# find the results that aren't nan
mask = ~np.isnan(result.Q)
result = result[mask]

# build the features
features = pd.DataFrame()
# mags first
features = result[['u','g','r', 'i', 'z']]

# colors
colors = combinations('zirgu', 2)
for i in colors:
    features['%s-%s' % (i[0], i[1])] = result.loc[:,i[0]] - result.loc[:,i[1]]

# colors squared
colors = combinations('zirgu', 2)
# for i in colors:
#     features['(%s-%s)2' % (i[0], i[1])] = (result.loc[:,i[0]] -
#             result.loc[:,i[1]])**2

# make datasets
X = features.values
y = result.Q.values
X_train, X_test, y_train, y_test = train_test_split(X,y,
        test_size=0.50)


# # use a full grid over all parameters -- for the random search
param_grid = {"max_depth": [3, None],
              "max_features": sp_randint(1, 15),
              "min_samples_split": sp_randint(1, 15),
              "min_samples_leaf": sp_randint(1, 10),
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"],
              "n_estimators": sp_randint(5, 100)}

scores = ['recall']

for score in scores:
    print("# Tuning hyper-parameters for %s" % score)
    print()

    # clf = GridSearchCV(RandomForestClassifier(),
    #    param_grid=param_grid,
    #    cv=5, scoring='%s_weighted' % score, n_jobs=-1)

    # perform a randomized grid search for the best possible 
    n_iter_search = 50
    clf = RandomizedSearchCV(ExtraTreesClassifier(),
            param_distributions=param_grid,
            n_iter=n_iter_search, cv=5, scoring='accuracy',
            n_jobs=-1)

    clf.fit(X, y)

    print("Best parameters set found on development set:")
    print()
    print(clf.best_params_)
    print()
    print("Grid scores on development set:")
    print()
    for params, mean_score, scores in clf.grid_scores_:
        print("%0.3f (+/-%0.03f) for %r"
             % (mean_score, scores.std() * 2, params))
    #print()

    # print("Detailed classification report:")
    # print()
    # print("The model is trained on the full development set.")
    # print("The scores are computed on the full evaluation set.")
    # print()
    # y_true, y_pred = y_test, clf.predict(X_test)
    # print(classification_report(y_true, y_pred))
    # print()
    # print(confusion_matrix(y_true, y_pred))

for j in range(10):
    print(j)
    magDict = {}
    with hdf.File('./truth/truth'+str(j).zfill(2)+'_Oii.hdf5', 'r') as f:
        dset = f['truth%s_Oii' % (str(j).zfill(2))]
        magDict['u'] = dset['OMAG'][:,0] # u band
        magDict['g'] = dset['OMAG'][:,1] # g band
        magDict['r'] = dset['OMAG'][:,2] # r band
        magDict['i'] = dset['OMAG'][:,3] # i band
        magDict['z'] = dset['OMAG'][:,4] # z band

    # we only want the g mag < 22 galaxies
    mask = np.where(magDict['g'] < 22)[0]
    print(mask.size)
    # create a data array for everything to fit into
    #features = -np.ones((mask.size, 25))
    features = -np.ones((mask.size, 15))

    # mags
    for i, m in enumerate('ugriz'):
        features[:,i] = magDict[m][mask]
    
    # colors
    colors = combinations('zirgu', 2)
    for i, c in enumerate(colors):
        features[:,i+5] = magDict[c[0]][mask] - magDict[c[1]][mask]
        # colors squared
        #features[:,i+15] = (magDict[c[0]][mask] - magDict[c[1]][mask])**2

    # now we make the predictions based on the new features we've created
    Qs = clf.predict(features)

    print(np.where(Qs == 0)[0].size/float(Qs.size))
    print(np.where(Qs == 1)[0].size/float(Qs.size))
    print(np.where(Qs == 2)[0].size/float(Qs.size))
        
    # with hdf.File('../truth/truth'+str(j).zfill(2)+'_Oii.hdf5', 'a') as f:
    #     values = -np.ones(magDict['u'].size)
    #     values[mask] = Qs
    #     f['Q'] = values


        