import pandas as pd
from glob import glob
import numpy as np
from astLib import astStats as ast
from astLib import astCalc as aca
from sklearn import mixture
import h5py as hdf
import emcee
import sys

# buzzard simulation cosmology
aca.H0 = 70
aca.OMEGA_M0 = 0.286
aca.OMEGA_L0 = 0.714

# millennium cosmology
#aca.H0 = 73
#aca.OMEGA_M0 = 0.25
#aca.OMEGA_L0 = 0.75


def findClusterCenterRedshift(data, errors=False):
    ''' Finds the center of the cluster in redshift space using the
    biweightlocation estimator. Puts the result in the ClusZ column of the data
    array. If errors = True, then this will return the 95% C.I. from 1000
    bootstrap shuffles.

    '''

    x = np.copy(data.redshift[data.interloper == 'NO'].values)
    w = np.copy(data.redshift_err.values) * 2.
    #avgz = ast.biweightLocation(x, tuningConstant=6.0)
    #return ast.biweightClipped(data['Z'], 6.0, 3)['biweightLocation']
    #return ast.biweightLocation(data['Z'], tuningConstant=6.0)
    avgz = np.average(x, weights=1. / w)

    #print len(x)
    if errors:
        ci = ast.bootstrap(x, np.average, weights=1. / w)
        return avgz, ci
    else:
        return avgz


def calcLOSVD(data, errors=False):
    ''' Using the previously computed LOSVs we will find the LOSVD. This will
    give us a few options to do that based on the number of objects that we
    have. If errors = True, then this will return the 95% C.I. from 1000
    bootstrap shuffles.

    '''

    if data.interloper.value_counts().NO >= 15:
        x = np.copy(data.LOSV[data.interloper == 'NO'].values)
        LOSVD = ast.biweightScale_test(x, tuningConstant=9.0)
        if errors:
            ci = ast.bootstrap(x, ast.biweightScale_test, tuningConstant=9.0)
            return LOSVD, ci
        else:
            return LOSVD

    else:
        x = np.copy(data.LOSV[data.interloper == 'NO'].values)
        LOSVD = ast.gapperEstimator(x)
        if errors:
            ci = ast.bootstrap(x, ast.gapperEstimator)
            return LOSVD, ci
        else:
            return LOSVD


def calc_mass_Evrard(data, A1D=1177., alpha=0.364):
    ''' This uses the relation from Munari2013 to calculate the halo mass from
    the observed velocity dispersion. The chosen scaling relations are from
    their table 1 which has been calibrated using galaxies and not dark matter
    halos only.

    '''

    avgz = findClusterCenterRedshift(data)
    vd = data.LOSVD.values[0]

    if avgz == None:
        pass
    else:
        return 1e15 / (aca.H0 * aca.Ez(avgz) / 100.) * (vd / A1D)**(1 / alpha)


def findLOSVDgmm(data):
    ''' Finds the LOSVD using a gaussian mixture model. Not currently used. '''

    LOSV = data['LOSV'][:, np.newaxis]
    lowest_bic = np.infty
    bic = []
    n_components_range = range(1, 4)
    cv_types = ['spherical', 'tied', 'diag', 'full']
    #cv_types = ['diag']
    for cv_type in cv_types:
        for n_components in n_components_range:
            # Fit a mixture of Gaussians with EM
            gmm = mixture.GMM(n_components=n_components,
                              covariance_type=cv_type,
                              n_init=10)
            gmm.fit(LOSV)
            bic.append(gmm.bic(LOSV))
            if bic[-1] < lowest_bic:
                lowest_bic = bic[-1]
                best_gmm = gmm

    # figure things out -- this comes from the wikipedia article on mixture
    # distributions. I have no idea if it is completely trustworthy.
    #covars = best_gmm.covars_.ravel()
    #weights = best_gmm.weights_.ravel()
    #means = best_gmm.means_.ravel()
    #wmeans = np.sum(weights*means)

    #parts = weights * ((means - wmeans)**2 + covars)
    #newLOSVD = np.sqrt(np.sum(parts))

    ## now we resample and then see
    dx = np.linspace(LOSV.min() - 100, LOSV.max() + 100, 1000)
    logprob, responsibilities = best_gmm.score_samples(dx[:, np.newaxis])
    pdf = np.exp(logprob)

    normedPDF = pdf / np.sum(pdf)

    u = np.sum(dx * normedPDF)
    data['LOSVDgmm'] = np.sqrt(np.sum(normedPDF * (dx - u)**2))
    return data


def findLOSVDmcmc(data):
    ''' Find the LOSVD and mean velocity using the MCMC and the likelihood
    function from walker2006. This tends to work better than any of the other
    methods and is the method we are using for the DES paper.

    '''

    def log_prior(theta, LOSV):
        sigma, mu = theta
        if not 50 < sigma < 1400:
            return -np.inf
        if not LOSV.min() < mu < LOSV.max():
            return -np.inf

        return 1

    def log_likelihood(theta, LOSV, LOSV_err):
        sigma, mu = theta
        #print(theta)
        # break long equation into three parts
        a = -0.5 * np.sum(np.log(LOSV_err**2 + sigma**2))
        b = -0.5 * np.sum((LOSV - mu)**2 / (LOSV_err**2 + sigma**2))
        c = -1. * LOSV.size / 2. * np.log(2 * np.pi)

        return a + b + c

    def log_posterior(theta, LOSV, LOSV_err):
        lp = log_prior(theta, LOSV)
        if not np.isfinite(lp):
            return -np.inf
        return lp + log_likelihood(theta, LOSV, LOSV_err)

    # get data
    LOSV = data['LOSV'].values
    try:
        LOSV_err = data['LOSV_err'].values
    except KeyError:
        LOSV_err = np.zeros_like(LOSV)

    ndim = 2  # number of parameters in the model
    nwalkers = 40  # number of MCMC walkers
    nburn = 50  # "burn-in" period to let chains stabilize
    nsteps = 300  # number of MCMC steps to take

    # set theta near the maximum likelihood, with
    np.random.seed()
    m = np.random.normal(np.mean(LOSV), scale=1, size=(nwalkers))
    s = np.random.normal(200, scale=1, size=(nwalkers))
    starting_guesses = np.vstack([s, m]).T

    sampler = emcee.EnsembleSampler(nwalkers,
                                    ndim,
                                    log_posterior,
                                    args=[LOSV, LOSV_err])
    sampler.run_mcmc(starting_guesses, nsteps)

    samples = sampler.chain[:, nburn:, :].reshape((-1, ndim))
    sigma_rec, mean_rec = map(
        lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
        zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    #data['LOSVD'] = sigma_rec[0]
    #data['LOSVD_err'] = sigma_rec[1], sigma_rec[2]

    #return data, samples

    return sigma_rec, samples


if __name__ == "__main__":
    ''' This takes the membership information from the member catalogs and
    computes the LOSVD and cluster mass from the member galaxies. Right now it
    is using a simple power law to do this, but things will change in the
    future when I add the stuff from the des study.

    '''

    if len(sys.argv) == 1:
        # make results container
        results = np.zeros(
            (10, ),
            dtype=[('ID', 'a', 25), ('SOURCES', '>i4'), ('Q0', '>i4'),
                   ('Q1', '>i4'), ('MEMBERS', '>i4'), ('Zc', '>f4'),
                   ('Zc_err', '>f4'), ('LOSVD', '>f4'), ('LOSVD_err', '>f4'),
                   ('LOSVD_dist', '>f4',
                    (10000, )), ('MASS', '>f4'), ('MASS_err', '>f4')])

        clusters = glob('./members/*_members.csv')
        results['ID'] = [c.rstrip('_members.csv').split('/')[-1]
                         for c in clusters]
        for i, c in enumerate(clusters):
            data = pd.read_csv(c)
            results['Q0'][i] = data[data.Q == 0].shape[0]
            results['Q1'][i] = data[data.Q == 1].shape[0]
            allobs = pd.read_csv('./redshifts/%s_redshifts.csv' %
                                 results['ID'][i])
            results['SOURCES'][i] = allobs[~np.isnan(allobs.Q)].shape[0]

            # filter out the non-members
            data = data[data.interloper == 'NO']
            results['MEMBERS'][i] = data.shape[0]

            # find redshift
            clusz = findClusterCenterRedshift(data, errors=True)
            results['Zc'][i] = clusz[0]

            # errors on a weighted mean, Eq. 4.19, Bevington
            results['Zc_err'][i] = 1. / np.sqrt(np.sum(
                1. / data.redshift_err.values))

            # calculate the LOSVD with the mcmc, without error to start
            LOSVD, LOSVD_dist = findLOSVDmcmc(data)
            results['LOSVD'][i] = LOSVD[0]
            results['LOSVD_dist'][i] = LOSVD_dist[:, 0]  # only LOSVD column
            results['LOSVD_err'][i] = np.std(LOSVD_dist[:, 0])

            # calculate the mass
            data['LOSVD'] = LOSVD[0]
            results['MASS'][i] = calc_mass_Evrard(data)

            # propigation of errors, Bevington
            results['MASS_err'][i] = results['MASS'][i]/0.364 *\
                        (results['LOSVD_err'][i]/ results['LOSVD'][i])

        results['MASS_err'] = 0.434 * results['MASS_err'] / results['MASS']
        results['MASS'] = np.log10(results['MASS'])

        with hdf.File('results_cluster.hdf5', 'w') as f:
            f['cluster props'] = results
            f.flush()

    elif len(sys.argv) == 2:
        cluster = './members/' + sys.argv[1] + '_members.csv'
        data = pd.read_csv(cluster)
        mass = calc_mass_Evrard(data)
        cluster = cluster.rstrip('_members.csv')
        print cluster, mass / 1e15
