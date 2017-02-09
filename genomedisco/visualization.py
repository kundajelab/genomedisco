
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

def plot_dds(dd_list,dd_names,out,approximation=10000):
    assert len(dd_list)==len(dd_names)
    
    fig, plots = plt.subplots(nrows=1, ncols=1)
    fig.set_size_inches(5, 5)
    colors=['red','blue']
    for dd_idx in range(len(dd_names)):
        dd_name=dd_names[dd_idx]
        dd=list(dd_list[dd_idx].values())
        x=list(dd_list[dd_idx].keys())
        sorted_x=np.argsort(np.array(x))
        x_plot=[]
        dd_plot=[]
        x_idx=0
        while x_idx<len(x):
            x_plot.append(x[sorted_x[x_idx]]*approximation)
            dd_plot.append(dd[sorted_x[x_idx]])
            x_idx+=1
        plots.plot(x_plot[1:],dd_plot[1:],c=colors[dd_idx],label=dd_names[dd_idx])
        plots.set_yscale('log',basey=10)
        plots.set_xscale('log',basex=10)
        plots.set_xlabel('Log10 Distance (bp)')
        plots.set_ylabel('Log10 Contact Probability')
    plots.legend(loc=3)
    fig.tight_layout()

    plt.savefig(out+'.png')
