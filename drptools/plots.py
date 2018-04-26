"""Some plotting utilities."""


from __future__ import print_function
import pylab
import IPython


def plot_coordinates(all_data, filtered_data, cluster_coord=None, title=None):
    """Plot a map of coordinates before and after a circular cut around the cluster center."""
    fig = pylab.figure()
    ax = fig.add_subplot(111, xlabel='ra', ylabel='dec')
    ax.scatter(all_data['coord_ra_deg'], all_data['coord_dec_deg'],
               color='k', label='All data', s=10)
    ax.scatter(filtered_data['coord_ra_deg'], filtered_data['coord_dec_deg'], color='b',
               label="Filtered data", s=10)
    if cluster_coord is not None:
        ax.scatter([cluster_coord[0]], [cluster_coord[1]], color='r', label='Cluster center',
                   marker='x', s=60)
    if title is not None:
        ax.set_title(title)
    ax.legend(loc='lower left', scatterpoints=1,
              frameon=False, fontsize='small')
    pylab.show()


def plot_patches(catalog, clust_coords=None):
    """Plot patches in the RA/DEC parameter space."""
    colors = pylab.cm.jet(pylab.linspace(0, 1, len(set(catalog['patch']))))
    fig = pylab.figure()
    ax = fig.add_subplot(111, xlabel="RA", ylabel="DEC")
    for i, path in enumerate(sorted(set(catalog['patch']))):
        filt = catalog['patch'] == path
        ra = catalog['coord_ra_deg'][filt]
        dec = catalog['coord_dec_deg'][filt]
        ax.scatter(ra, dec, color=colors[i], label=path, s=10)
        ax.vlines(min(ra), min(dec), max(dec), color='k')
        ax.vlines(max(ra), min(dec), max(dec), color='k')
        ax.hlines(min(dec), min(ra), max(ra), color='k')
        ax.hlines(max(dec), min(ra), max(ra), color='k')
        ax.legend(loc='best', numpoints=1, frameon=False)
    if clust_coords is not None:
        ax.scatter(clust_coords['ra'], clust_coords['dec'],
                   s=100, marker='s', color='k')
    pylab.show()


def js9(myimage):
    html = """
    <head>
    <link type="text/css" rel="stylesheet" href="https://js9.si.edu/js9/js9support.css">
    <link type="text/css" rel="stylesheet" href="https://js9.si.edu/js9/js9.css">
    <script type="text/javascript" src="https://js9.si.edu/js9/js9prefs.js"></script>
    <script type="text/javascript" src="https://js9.si.edu/js9/js9support.min.js"></script>
    <script type="text/javascript" src="https://js9.si.edu/js9/js9.js"></script>
    <script type="text/javascript" src="https://js9.si.edu/js9/js9plugins.js"></script>
    <link type="text/css" rel="stylesheet" href="https://js9.si.edu/js9/js9-allinone.css">
    <script type="text/javascript" src="https://js9.si.edu/js9/js9-allinone.js"></script>
    </head>
    <body onload="JS9.Preload('IMAGETOLOAD', {scale: 'log', zoom: 'to fit'});">
    <div class="JS9Menubar" data-width="600px"></div>
    <div class="JS9Toolbar" data-width="600px"></div>
    <div class="JS9" data-width="600px" data-height="600px"></div>
    <div style="margin-top: 2px;" data-width="600px">
    <div class="JS9Colorbar" data-width="600px"></div>
    href='javascript:JS9.Load("fits/casa.fits");'
    </div>
    </body>
    """.replace('IMAGETOLOAD', myimage)
    return IPython.display.HTML(html)