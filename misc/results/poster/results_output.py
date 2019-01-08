import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def print_table(df,display_table=True):
    int_str = '{:d}'.format
    sci_str = '{:.2e}'.format
    float_str = '{:.4f}'.format
    fstr = {'p' : int_str,
        'iterations' : int_str,
        'residual' : sci_str,
        'error' : sci_str,
        'domain_init' : float_str,
        'matrix_setup' : float_str,
        'amgx_setup' : float_str,
        'linsys_setup' : float_str,
        'linsolve' : float_str,
        'patch_solve' : float_str,
        'complete_solve' : float_str,
        'complete_setup' : float_str}

    styler = df.style
    styler = styler.set_properties(subset=None,width='100px')
    styler = styler.set_properties(**{'subset':'p','width': '50px', 'text-align': 'right'})
    styler = styler.set_properties(**{'subset':['residual','error'],'width': '80px', 'text-align': 'right'})
    styler = styler.set_properties(subset=['p','iterations'],width='30px')
    styler = styler.format(fstr)
    if display_table:
        display(styler)
    return styler
    
def fix_xticks(procs):
    
    p0 = np.log2(procs[0])
    p1 = np.log2(procs[-1])
    plt.xlim([2**(p0-1), 2**(p1+1)])
    
    pstr = ([str(p) for p in procs])
    plt.xticks(procs,pstr)
    # Suppress minor tick marks;  doesn't make sense for processor counts
    ax = plt.gca()
    ax.xaxis.set_minor_locator(plt.NullLocator())

def set_xticks(tick_major,tick_minor,fmt=None,units='seconds'):
    if units == 'seconds':
        scale_factor = 1
    elif units == 'minutes':
        scale_factor = 60
    elif units == 'hours':
        scale_factor = 3600

    if fmt is None:
        fmt_func = lambda value,tick_number : '{:.1f}'.format(value/scale_factor)
    else:
        fmt_func = lambda value,tick_number : fmt(value/scale_factor,tick_number)

    ax = plt.gca()
    ax.xaxis.set_major_locator(plt.MultipleLocator(scale_factor*tick_major))   # Major tick marks every 10 seconds
    ax.xaxis.set_minor_locator(plt.MultipleLocator(scale_factor*tick_minor))    # Minor tick marks every 5 seconds
    ax.xaxis.set_major_formatter(plt.FuncFormatter(fmt_func))  # Format tick marks

    plt.xlabel('Time ({:s})'.format(units))


    
def strong_scaling(df,field='walltime'):
    df.plot(x='p',y=field,logx=True,logy=True,style='.-',markersize=15,label=field.capitalize())
    procs = df['p'].values

    # Plot best-fit speed-up line
    t_strong = np.array(df[field].values)
    c = np.polyfit(np.log(procs[:-2]),np.log(t_strong[:-2]),1)
    plt.loglog(procs,np.exp(np.polyval(c,np.log(procs))),'r-',label='Best-fit (slope={:6.2f})'.format(c[0]),linewidth=1)
    c[0] = -1
    plt.loglog(procs,np.exp(np.polyval(c,np.log(procs))),'k--',label='Theoretical',linewidth=0.5)

    fix_xticks(procs)

    plt.legend()
    plt.show()
    

def efficiency(df,field='walltime'):
 
    # Efficiency
    plt.figure()
    procs = df['p'].values
    T0 = df[field][0]
    S = T0/df[field]
    E = 100*S/procs

    plt.semilogx(procs,E,'.-',markersize=15)
    plt.semilogx(procs,[100]*len(df),'k--',linewidth=2)

    plt.xlabel('p',fontsize=16)
    plt.ylabel('Efficiency (%)',fontsize=16)
    # plt.title("Efficiency (%)");
    plt.legend([field, 'Perfect efficiency'])
    
    fix_xticks(procs)
    
    plt.ylim([-1,110])
    plt.grid()
    plt.show()

    
def bar_plot(df,cols=None):
    if cols is None:
        cols = df.columns.values
    
    # For plotting
    df_plot = df[cols].iloc[::-1].copy()

    c = []
    for p in df['p'].values:
        c += ['{:d} proc(s)'.format(p)]
        
    df_plot.index = reversed(c)

    ax = df_plot.plot.barh(width=0.75)
    plt.xlabel('Time (seconds)');

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(reversed(handles), reversed(labels), loc='lower right')

    fmt = lambda value,tick_number : '{:.1f}'.format(value)
    set_xticks(10,5,fmt)                                 
    plt.grid()
    # plt.xlim([0,135*60]);
    plt.show()    

    
def barh_plot(df,cols=None):
    # For plotting (iloc[::-1] reverses order of the rows)
    if cols is None:
        cols = df.columns.values
    
    df_plot = df[cols].iloc[::-1].copy()
    c = []
    for p in df['p'].values:
        c += ['{:d} proc(s)'.format(p)]

    df_plot.index = reversed(c)
    
    ax = df_plot.plot.barh(width=0.85,stacked=True)
    plt.xlabel('Time (seconds)');

    handles, labels = ax.get_legend_handles_labels()
    # ax.legend(reversed(handles), reversed(labels), loc='lower right')

    fmt = lambda value,tick_number : '{:.1f}'.format(value)
    set_xticks(10,5,fmt)

    plt.grid()
    # plt.xlim([0,135*60]);
    plt.show()        