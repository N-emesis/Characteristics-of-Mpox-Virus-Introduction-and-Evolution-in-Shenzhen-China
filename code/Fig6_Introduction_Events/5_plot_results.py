import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import warnings

def set_style():
    # Defines colors based on the reference script
    global greys
    greys = ["#2E3440", "#3b4252", "#434C5E", "#4C566A"]
    global mpl
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica', 'Liberation Sans', 'sans-serif']
    mpl.rcParams['font.weight'] = 'light'
    mpl.rcParams['text.color'] = greys[0]
    mpl.rcParams['axes.labelcolor'] = greys[0]
    mpl.rcParams['xtick.color'] = greys[0]
    mpl.rcParams['ytick.color'] = greys[0]
    mpl.rcParams['figure.titlesize'] = 16
    mpl.rcParams['axes.titlesize'] = 16
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16
    mpl.rcParams['axes.edgecolor'] = greys[0]

    global sns
    if hasattr(sns.categorical, '_Old_Violin'):
        pass
    else:
        try:
            sns.categorical._Old_Violin = sns.categorical._ViolinPlotter
            class _My_ViolinPlotter(sns.categorical._Old_Violin):
                def __init__(self, *args, **kwargs):
                    super(_My_ViolinPlotter, self).__init__(*args, **kwargs)
                    self.gray = greys[0]
            sns.categorical._ViolinPlotter = _My_ViolinPlotter
        except AttributeError:
            pass

def plot_rarefaction_results(rarefaction_results, outname):
    # Filter only Shenzhen data
    shenzhen_results = rarefaction_results[rarefaction_results['region'] == 'Shenzhen'].copy()
    
    # +1 logic from plot.py reference
    shenzhen_results['importations'] = shenzhen_results['importations'] + 1

    inner_colors = {'Shenzhen': '#6483a4'}
    
    fig, axs = plt.subplots(1, 1, figsize=(12, 6), constrained_layout=True)
    
    sns.violinplot(x='n', y='importations', hue='region', cut=0, 
                   data=shenzhen_results, ax=axs,
                   palette=inner_colors, saturation=1, scale='width', width=0.7)
    
    # Customizing the inner boxplot drawn by seaborn
    # We want to make sure the median is a distinct white dot with a black border
    for line in axs.lines:
        if line.get_color() in ['w', 'white', '#ffffff'] or line.get_markeredgecolor() in ['w', 'white', '#ffffff'] or line.get_markerfacecolor() == 'black':
            line.set_marker('o')  # Ensure it is a circle
            line.set_color('white') # Set line color
            line.set_markerfacecolor('white') # Set face to white
            line.set_markeredgecolor('black') # Set edge to black
            line.set_markersize(8) # Size of the dot
            line.set_markeredgewidth(1.5) # Border width
    for c in axs.collections:
        if hasattr(c, 'get_facecolor') and len(c.get_facecolor()) > 0:
            if (c.get_facecolor()[0][:3] == [1., 1., 1.]).all() or (c.get_facecolor()[0][:3] == [0., 0., 0.]).all():
                c.set_facecolor('white') # Fallback for some seaborn versions
                c.set_edgecolor('black')
                c.set_linewidth(1.5)
                c.set_sizes([64]) 
                
    if axs.get_legend() is not None:
        axs.get_legend().remove()
        
    axs.set_xlabel('Number of samples', fontname='Arial', fontsize=24)
    axs.set_ylabel('Estimated number of\nintroduction events', fontname='Arial', fontsize=24)
    
    # Set tick labels font
    axs.set_xticklabels(axs.get_xticklabels(), fontname='Arial', fontsize=22)
    axs.set_yticklabels(axs.get_yticklabels(), fontname='Arial', fontsize=22)

    fig.savefig(f'{outname}.pdf')
    fig.savefig(f'{outname}.svg')
    plt.close(fig)

def main():
    import argparse
    set_style()
    parser = argparse.ArgumentParser()
    parser.add_argument('--rarefaction_dir', default='results/rarefaction/*/*_rarefaction.csv')
    args = parser.parse_args()
    
    rarefaction_results = pd.DataFrame()
    for file in glob.glob(args.rarefaction_dir):
        df = pd.read_csv(file)
        if isinstance(df['region'].iloc[0], list):
            df['region'] = df['region'].apply(lambda x: ''.join(x) if isinstance(x, list) else x)
        rarefaction_results = pd.concat([rarefaction_results, df], ignore_index=True)
        
    if len(rarefaction_results) == 0:
        print("No rarefaction results found.")
        return
        
    print(rarefaction_results['region'].unique())
    
    # Ensure n is numeric and sorted, limit to 180
    rarefaction_results['n'] = pd.to_numeric(rarefaction_results['n'])
    rarefaction_results = rarefaction_results[rarefaction_results['n'] <= 197]
    rarefaction_results = rarefaction_results.sort_values(by=['n', 'region'])
    
    os.makedirs('figures', exist_ok=True)
    
    plot_rarefaction_results(rarefaction_results, 'figures/rarefaction1_Shenzhen_Only')
    print(f"Successfully generated Shenzhen-only rarefaction plots: figures/rarefaction1_Shenzhen_Only.pdf and .svg")

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()