# -----------------------------------------------------------------------------
# Change input parameters below, and run the whole file in your favorite
# Python editor. If you set wd to 'current', the weight files will be saved
# to your current working directory.
# -----------------------------------------------------------------------------

# Input parameters
resolution = 120  # Change to the resolution of the input raster for r.mfilter
wd = "current"  # Where do you want to save your output files?
view_plot = True  # Do you want to view the kernel functions?
radius_QHF = 450  # Only change if you have good reasons to do so
radius_QHS = 210  # Only change if you have good reasons to do so

# -----------------------------------------------------------------------------
# Do not change things below, unless you know what you do
# -----------------------------------------------------------------------------


def createFilterFile(
    resolution,
    radius,
    file,
    title,
    exponent=1,
    wd="current",
    view_plot=True,
    save_plot=False,
):
    """
    Returns a file with inverse distance weight filter that can be
    used as input for the GRASS GIS function r.mfilter. In addition,
    it prints a plot to screen showing the shape of the inverse
    distance weight function.
    Todo: circular neighborhood option.

    Parameters
    ----------
    resolution : integer
        resolution to be used (should be the resolution of the current region).
    radius : integer
        Radius of neighbourhood in meters
    file : TYPE
        Name of output file
    title : string
        A one-line TITLE for the filter, used to construct a TITLE for the
        resulting raster map layer.
    exponent : TYPE
        Distance weighting exponent. The Default =1, which gives a inverse
        linear weight function. An exponent of 2 will give an inverse distance
        squared weighting
    wd : TYPE, optional
        Path to working directory. The default is 'current'.
    save_plot : Boolean (True/False)
       Save the plot to file. The default is False.

    Returns
    -------
    None.

    """

    # Libraries
    import os
    import math
    import numpy as np
    import matplotlib.pyplot as plt

    if wd != "current":
        try:
            os.chdir(wd)
            print("Working directory set to:  newPath")
        except OSError:
            print("Directory ", wd, "does not exist.")

    # Create weight-distance array

    radcel = math.ceil(radius / resolution)
    numcel = radcel * 2 + 1
    bn1 = (np.array(range(1, radcel + 1)) / (radcel + 1)) ** exponent
    bn = np.concatenate((bn1, np.array([1.0]), np.flip(bn1))).round(decimals=10)

    # Create matrix
    mfilter = np.zeros([numcel, numcel])
    for i in range(0, numcel):
        mfilter[i : (numcel + 1 - i - 1), i] = bn[i]
        mfilter[i : (numcel - i), numcel - (i + 1)] = bn[i]
        mfilter[i, i : (numcel + 1 - i - 1)] = bn[i]
        mfilter[numcel - (i + 1), i : (numcel - i)] = bn[i]

    # Header and footer
    header = "TITLE {}\nMATRIX {}".format(title, numcel)
    footer = "DIVISOR 0\nTYPE P"

    # Print to file
    np.savetxt(
        file,
        mfilter,
        fmt="%0.6g",
        header=header,
        comments="",
        footer=footer,
        delimiter=" ",
    )

    # Create figure
    if view_plot:
        steps = list(range(-radcel, 0)) + [0] + list(range(1, radcel + 1))
        distance = [i * resolution for i in steps]
        fig = plt.figure(1, figsize=(6, 4))
        plt.plot(distance, bn)
        plt.title("inverse distance weight function")
        plt.xlabel("distance (meter)")
        plt.ylabel("weight")
        filename = os.path.splitext(file)[0] + ".png"
        if save_plot:
            fig.savefig(filename, dpi=300)
        plt.show()


exponent = 1
createFilterFile(
    resolution=resolution,
    radius=radius_QHF,
    file=f"QHFwaf_{resolution}.txt",
    title="distance weighted average QHF",
    exponent=exponent,
    view_plot=view_plot,
    wd=wd,
)
createFilterFile(
    resolution=resolution,
    radius=radius_QHS,
    file=f"QHSwaf_{resolution}.txt",
    title="distance weighted average QHS",
    exponent=exponent,
    view_plot=view_plot,
    wd=wd,
)
