import pandas as pd
from astLib import astCoords
from numpy import round, isnan
import sys

def hms(ra):
    try:
        return astCoords.decimal2hms(ra, ':')
    except UnboundLocalError:
        return 'NaN'

def dms(dec):
    try:
        return astCoords.decimal2dms(dec, ':')
    except UnboundLocalError:
        return 'NaN'

def main(cluster):
    df = pd.read_csv('./../members/'+cluster+'_members.csv')
    df_part = df[['tile', 'dither', 'fiber', 'ra', 'dec', 'r',
    'redshift','redshift_err', 'Q', 'interloper', 'separation', 'LOSV',
    'LOSV_err']]

    # remove the un-observed galaxies
    x = ~isnan(df_part['redshift'])
    df_part = df_part[x]

    # now we are going to round some of the columns to the right length
    df_part[['r']] = df_part[['r']].apply(round, args=(2,))
    df_part[['redshift']] = df_part[['redshift']].apply(round, args=(4,))
    df_part[['separation']] = df_part[['separation']].apply(round, args=(2,))

    # now we have to convert the ra/dec to sexigesimal
    df_part[['ra']] = [hms(ra) for ra in df_part[['ra']].values.ravel()]
    df_part[['dec']] = [dms(dec) for dec in df_part[['dec']].values.ravel()]

    # replace the interloper with checkmarks
    df_part.interloper[df_part.interloper == 'NO'] = r'$\checkmark$'
    df_part.interloper[df_part.interloper == 'YES'] = r'...'


    df_part.to_latex(cluster+'_table.tex', index=False)


if __name__ == "__main__":
    main(sys.argv[1])

