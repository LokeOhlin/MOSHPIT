def fmt_10(x, base):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


color_scheme = {
    "dust" : "#E69F00",
    "dust_analytic" : "#F0E442",
    "gas"  : "#0072B2",
    "gas_analytic": "#56B4E9"
        }

label_scheme = {
    "rhodust" : r"\rho_{d}", 
    "vdust" : r"v_{r}",

    "rhogas" : r"\rho_{g}", 
    "vgas" : r"u_{r}",
    "eint" : r"\epsilon",
    "etot" : r"e",
    "pres" : r"p"

}

one_width = 3.34

xlabel_size = one_width*0.175
ylabel_size = one_width*0.175

ax_width  = one_width - 2*ylabel_size
ax_height = 4*ax_width/5

min_sep = 0.00



def get_figure_parameters(ncols, nrows, sharex = False, sharey = False, ax_height_ratio = 4/5):
    ax_height = ax_width* ax_height_ratio
    one_height = ax_height + 2*xlabel_size
    
    if sharex :
        fig_height = nrows*(ax_height + min_sep) + 2*xlabel_size - min_sep
    else:
        fig_height = nrows*(ax_height + xlabel_size + 0.5*min_sep) + xlabel_size - 0.5*min_sep

    if sharey :
        fig_width = ncols*(ax_width + min_sep) + 2*ylabel_size - min_sep
    else:
        fig_width = ncols*(ax_width + ylabel_size + 0.5*min_sep) + ylabel_size - 0.5*min_sep


    bottom = xlabel_size/fig_height
    top = 1-bottom
    left  = ylabel_size/fig_width
    right = 1-left

    if sharex :
        hspace = min_sep/ax_height
    else:
        hspace = (xlabel_size+0.5*min_sep)/ax_height

    if sharey :
        wspace = min_sep/ax_width
    else:
        wspace = (ylabel_size+0.5*min_sep)/ax_width

    pars = {"top": top,
            "bottom": bottom,
            "left" : left,
            "right" : right,
            "hspace":hspace,
            "wspace":wspace}


    return (fig_width, fig_height), pars


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    figsize, pars = get_figure_parameters(1,1)
    fig, ax = plt.subplots(1,1, figsize = figsize)
    plt.subplots_adjust(**pars)
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")

    figsize, pars = get_figure_parameters(2,2, sharex = True, sharey = True)
    fig, axes = plt.subplots(2,2, figsize = figsize, sharex = True, sharey = True)
    plt.subplots_adjust(**pars)
    for ax in axes.ravel():
        ax.set_xlim(-1,1)
        ax.set_ylim(-1,1)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")

    plt.show()

