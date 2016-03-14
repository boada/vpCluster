import pandas as pd
from glob import glob
import numpy as np
from itertools import combinations
from sklearn.svm import SVC
from sklearn.cross_validation import train_test_split
from sklearn.grid_search import GridSearchCV, RandomizedSearchCV
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from scipy.stats import randint as sp_randint


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
for i in colors:
    features['(%s-%s)2' % (i[0], i[1])] = (result.loc[:,i[0]] -
            result.loc[:,i[1]])**2

# make datasets
X = features.values
y = result.Q.values
X_train, X_test, y_train, y_test = train_test_split(X,y,
        test_size=0.5)

# Set the parameters by cross-validation
param_grid = [{'kernel': ['rbf'],
                    'gamma': [1e-3, 1e-4],
                     'C': [1, 10, 100, 1000]},
                    {'kernel': ['linear'],
                    'C': [1, 10, 100, 1000]}]

scores = ['precision', 'recall']

for score in scores:
    print("# Tuning hyper-parameters for %s" % score)
    print()

    clf = GridSearchCV(SVC(), param_grid=param_grid, cv=5,
                       scoring='%s_weighted' % score)
    clf.fit(X_train, y_train)

    print("Best parameters set found on development set:")
    print()
    print(clf.best_params_)
    print()
    print("Grid scores on development set:")
    print()
    for params, mean_score, scores in clf.grid_scores_:
        print("%0.3f (+/-%0.03f) for %r"
            % (mean_score, scores.std() * 2, params))
    print()

    print("Detailed classification report:")
    print()
    print("The model is trained on the full development set.")
    print("The scores are computed on the full evaluation set.")
    print()
    y_true, y_pred = y_test, clf.predict(X_test)
    print(classification_report(y_true, y_pred))
    print()
    print(confusion_matrix(y_true, y_pred))


# # use a full grid over all parameters -- for the random search
param_grid = {"max_depth": [3, None],
              "max_features": sp_randint(1, 25),
              "min_samples_split": sp_randint(1, 25),
              "min_samples_leaf": sp_randint(1, 10),
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"],
              "n_estimators": sp_randint(5, 100)}

# use a full grid over all parameters -- for the grid search
# param_grid = {"max_depth": [3, None],
#                 "max_features": [1, 3, 10, 15, 25],
#                 "min_samples_split": [1, 3, 10, 15, 25],
#                 "min_samples_leaf": [1, 3, 10],
#                 "bootstrap": [True, False],
#                 "criterion": ["gini", "entropy"],
#                 "n_estimators": [5, 20, 50, 75]}

scores = ['precision', 'recall']

for score in scores:
    print("# Tuning hyper-parameters for %s" % score)
    print()

    # clf = GridSearchCV(RandomForestClassifier(),
    #    param_grid, cv=5, scoring='%s_weighted' % score)

    n_iter_search = 50
    clf = RandomizedSearchCV(RandomForestClassifier(),
            param_distributions=param_grid,
            n_iter=n_iter_search, cv=5, scoring='%s_weighted' % score,
            n_jobs=-1)

    clf.fit(X_train, y_train)

    print("Best parameters set found on development set:")
    print()
    print(clf.best_params_)
    print()
    print("Grid scores on development set:")
    print()
    #for params, mean_score, scores in clf.grid_scores_:
        #print("%0.3f (+/-%0.03f) for %r"
        #      % (mean_score, scores.std() * 2, params))
    #print()

    print("Detailed classification report:")
    print()
    print("The model is trained on the full development set.")
    print("The scores are computed on the full evaluation set.")
    print()
    y_true, y_pred = y_test, clf.predict(X_test)
    print(classification_report(y_true, y_pred))
    print()
    print(confusion_matrix(y_true, y_pred))



### This is the best params found by the random search.. they seem better than
### the grid search.

params = {'bootstrap': True, 'min_samples_leaf': 3, 'n_estimators': 91,
'min_samples_split': 6, 'criterion': 'entropy', 'max_features': 12,
'max_depth': 3}



########################
# this stuff isn't used.
#######################
# # Create the RFE object and compute a cross-validated score.
# svc = SVC(kernel="rbf")
# # The "accuracy" scoring is proportional to the number of correct
# # classifications
# rfecv = RFECV(estimator=svc, step=1, cv=StratifiedKFold(y_train, 2),
#                       scoring='accuracy')
# rfecv.fit(X_train, y_train)
#
# print("Optimal number of features : %d" % rfecv.n_features_)


