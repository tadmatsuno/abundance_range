import solara
import matplotlib.pyplot as plt
import matplotlib
import pandas
import numpy as np
from scipy import stats
from astropy.table import Table


matplotlib.rcParams['font.size'] = 12
all_element = ['Li','C','N','O','Na','Mg',\
               'Al','Si','Ca','Ti','V',
               'Cr','Mn','Ni','Co','Sr',
               'Y','Zr','Ba','Eu']

sagamp = pandas.read_csv('./4MOST_investigate_abundancerange0209MP.tsv',sep='\t')
sagamp_additional = {'Na':'./4MOST_investigate_abundancerange0209MP_Na.tsv'}
sagamr = pandas.read_csv('./4MOST_investigate_abundancerange0209MR.tsv',sep='\t')
sagamr_additional = {'Na':'./4MOST_investigate_abundancerange0209MR_Na.tsv'}
#galah = Table.read('./GALAHDR3_snr50_spfeh_flags.fits')


quantiles_print = [5,95]

def plot_one_axis(elem,ax,MP=False,ranges=None,SAGA=False,GALAH=False,APOGEE=False):
    if SAGA:
        if MP:
            if elem in sagamp_additional.keys():
                data = pandas.read_csv(sagamp_additional[elem],sep='\t')
            else:
                data = sagamp
        else:
            if elem in sagamr_additional.keys():
                data = pandas.read_csv(sagamr_additional[elem],sep='\t')
            else:
                data = sagamr
        col = f'[{elem}/Fe]'
        if 'Li' in elem:
            col = f'log-e(Li)'
        val = data[col]
        if val.dtypes == 'object':
            mask = ~(val.str.startswith('<').astype(bool)|\
                            val.str.startswith('>').astype(bool)|\
                            val.str.startswith('.').astype(bool)|\
                            val.str.contains('[a-zA-Z]',regex=True).astype(bool))
            dsmall = data[mask]
        else:
            dsmall = data.copy()
        mask2 = np.isfinite(dsmall[col].astype(float))
        dsmall = dsmall[mask2]
        dsmall = dsmall[~dsmall['# Object'].duplicated(keep='first')]
        val_num = dsmall[col].astype(float)
        if dsmall['[Fe/H]'].dtype == 'object':
            feh = np.where(dsmall['[Fe/H]'].str.startswith('<'),np.nan,dsmall['[Fe/H]']).astype(float)
        else:
            feh = dsmall['[Fe/H]']
        ax.plot(feh,val_num,'ko',ms=1.,alpha=0.5)
        if ranges:
            xmin,xmax = ax.get_xlim()
            lpercentile = stats.percentileofscore(val_num,ranges[elem][0])
            upercentile = stats.percentileofscore(val_num,ranges[elem][1])
            ax.text(xmin+0.05*(xmax-xmin),ranges[elem][0]-0.05*(ranges[elem][1]-ranges[elem][0]),
                    f'{lpercentile:.1f}',color='k',verticalalignment='top',horizontalalignment='left')
            ax.text(xmin+0.05*(xmax-xmin),ranges[elem][1]+0.05*(ranges[elem][1]-ranges[elem][0]),
                    f'{upercentile:.1f}',color='k',verticalalignment='bottom',horizontalalignment='left')
        
                    
@solara.component
def Page():
    elements,set_elements = solara.use_state(all_element)
    mpmr,set_mpmr = solara.use_state('MR')
    ranges = {}
    set_ranges = {}
    for elem in all_element:
        if elem == 'Li':
            ranges[elem],set_ranges[elem] = solara.use_state([0,3])
        else:
            ranges[elem],set_ranges[elem] = solara.use_state([-1,1])


    nelem = len(elements)
    if nelem <= 2:
        nrow = 1
        ncol = np.maximum(nelem,1)
    elif nelem <=4:
        nrow = 2
        ncol = 2
    elif nelem <= 6:
        nrow = 2
        ncol = 3
    elif nelem <= 9:
        nrow = 3
        ncol = 3
    elif nelem <= 12:
        nrow = 3
        ncol = 4
    elif nelem <= 16:
        nrow = 4
        ncol = 4
    elif nelem <= 20:
        nrow = 4
        ncol = 5

    fig, axs = plt.subplots(nrow, ncol, 
        figsize=(12, 12),sharex=True)
    axs = np.atleast_1d(axs).ravel()
    
    for ii,elem in enumerate(elements):
        ax = axs[ii]
        ax.set_title(elem)
        if mpmr == 'MP':
            ax.set_xlim(-5,-2)
        else:
            ax.set_xlim(-2,0.5)
        ax.set_ylim(-1.5,1.5)
        if 'Li' in elem:
            ax.set_ylim(-0.5,5)
        elif elem in ['C','N','O']:
            ax.set_ylim(-1.5,4)
        elif elem in ['Mg','Si','Ca','Ti']:
            ax.set_ylim(-1.,2.5)
        elif elem in ['Sr','Y','Zr','Ba','Ce','Eu']:
            ax.set_ylim(-3.,4)
        plot_one_axis(elem,ax,MP=mpmr=='MP',ranges=ranges,SAGA=True)
        xmin,xmax = ax.set_xlim()
        ax.fill_between([xmin,xmax],
            [ranges[elem][0]]*2,[ranges[elem][1]]*2,color='C7',alpha=0.3)

    with solara.Column() as main:
        with solara.Columns([1,2]):
            with solara.Card():
                with solara.Row():
                    solara.ToggleButtonsSingle(value=mpmr, values=['MP','MR'], on_value=set_mpmr)
                solara.Text('\nSelection of elements')
                solara.SelectMultiple('Elements',
                        values=elements,
                        all_values=all_element,
                        on_value=set_elements)
#                solara.ToggleButtonsMultiple(value=elements,
#                        values=all_element,
#                        on_value=set_elements,dense=True)
                with solara.Row():
                    solara.Button('Select all',on_click=lambda: set_elements(all_element))
                    solara.Button('Clear all',on_click=lambda: set_elements([]))
                solara.Text('\nRanges for individual elements')
                for elem in all_element:
                    solara.SliderRangeFloat(f'{elem:2s} range',min=-4,max=4,step=0.1,
                        value=ranges[elem],on_value= set_ranges[elem],thumb_label='always')

            fig.tight_layout()
            solara.FigureMatplotlib(fig,dependencies=[ranges,mpmr])
#            solara.FigureMatplotlib(fig,dependencies=[elements,ranges,mpmr])
Page()
