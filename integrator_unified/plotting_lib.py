import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.ticker import LogLocator, FuncFormatter, MultipleLocator
from IPython.display import display

import matplotlib.pyplot as plt
from analysis_lib import *

# default colormap used by the plotting functions (can be overridden elsewhere)
cmap = mpl.cm.get_cmap("tab10")


def geodesic_plot(params, show=False, save=False, xaxis=None, range=[0, -1]):
    filename = give_filename(params, tidal=False)
    data = load_geodesic(filename)

    spin = params["spin"]
    charge = params["charge"]
    x = params["x"]

    if xaxis:
        xdata = data[xaxis]
    else:
        xdata = data.index

    # Plot r in log scale
    plt.figure(figsize=(9, 6))
    # plt.plot(xdata, data['r'], color=cmap.colors[0])
    plt.plot(
        xdata[range[0] : range[1]], data["r"][range[0] : range[1]], color=cmap.colors[0]
    )
    plt.xlabel("Iteration")
    plt.ylabel("r")
    plt.yscale("log")

    # Set background colour to white
    plt.gca().patch.set_facecolor("white")

    # Draw r=x line
    plt.axhline(y=x, color=cmap.colors[1], linestyle="--")
    # Draw horizon if spin^2 + charge^2 < 1
    if spin**2 + charge**2 < 1:
        horizon = 1 + np.sqrt(1 - spin**2 - charge**2)
        plt.axhline(y=horizon, color=cmap.colors[2], linestyle="--")
    else:
        horizon = None
    # Draw repelling line if charge != 0
    if charge != 0:
        rep = charge * charge / 2
        plt.axhline(y=rep, color=cmap.colors[3], linestyle="--")

    # Add legend at top right
    legend_list = ["r", "IBSO = {:.2f}".format(x)]
    if charge != 0:
        legend_list.append("Inner root = {:.2f}".format(rep))
    if horizon:
        legend_list.append("Horizon = {:.2f}".format(horizon))

    plt.legend(legend_list, loc="best")

    plt.tight_layout()

    if save:
        plt.savefig(give_filename(params, tidal=False).replace(".dat", ".png"), dpi=500)

    if show:
        plt.show()


def tidal_plot(
    params,
    show=False,
    save=False,
    legend=True,
    dpi=1000,
    figsize=(10, 7),
    fontsize=18,
    cutoff=None,
    major_length=3,
    minor_length=2,
    xaxis=None,
    custom_cmap=None,
):
    filename = give_filename(params, tidal=False)
    data = load_geodesic(filename, cutoff=cutoff)

    filename = give_filename(params, tidal=True)
    tidal_data = load_tensor_cpp(filename, cutoff=cutoff)

    index_i = 0
    if cutoff:
        index_f = cutoff
    else:
        try:
            index_f = stability[1] - 50
        except:
            index_f = len(data)

    # Plot r, theta, and second column of tidal tensor on the same plot
    fig, ax1 = plt.subplots(figsize=figsize)

    # Set background colour to white
    fig.patch.set_facecolor("white")

    ax1.plot(
        data.index[index_i:index_f],
        data["r"][index_i:index_f],
        color=cmap.colors[0],
        linestyle="--",
    )
    ax1.set_ylim(0.3, 1e6 * 3)
    ax1.set_yscale("log")  # Fix ax1 y-axis limits

    ax1.set_xlabel("Iteration", fontsize=fontsize)
    ax1.set_ylabel("$r$", color=cmap.colors[0], fontsize=fontsize)
    ax1.spines["left"].set_color(cmap.colors[0])
    ax1.tick_params(
        "y",
        which="major",
        colors=cmap.colors[0],
        length=major_length,
        width=1,
        labelsize=fontsize,
    )
    ax1.tick_params(
        "y",
        which="minor",
        colors=cmap.colors[0],
        length=minor_length,
        width=1,
        labelsize=fontsize,
    )
    ax1.tick_params(
        "x", which="major", length=major_length, width=1, labelsize=fontsize
    )
    ax1.tick_params(
        "x", which="minor", length=minor_length, width=1, labelsize=fontsize
    )
    # Only tick even powers of 10
    ax1.yaxis.set_major_locator(LogLocator(base=100))
    ax1.minorticks_off()
    # Put x-label at the top
    # ax1.xaxis.set_label_position("top")
    # Customise y-axis label position
    ax1.yaxis.set_label_coords(-0.02, 0.5)
    # Customise x-axis label position
    # ax1.xaxis.set_label_coords(0.5, -0.05)

    ax2 = ax1.twinx()

    ax2.plot(
        data.index[index_i:index_f],
        data["theta"][index_i:index_f],
        color=cmap.colors[1],
        linestyle=":",
    )
    ax2.set_ylim(0, np.pi)
    ax2.set_ylabel("$\\theta$", color=cmap.colors[1], fontsize=fontsize)
    ax2.spines["right"].set_color(cmap.colors[1])
    ax2.spines["left"].set_color(cmap.colors[0])
    ax2.tick_params(
        "y",
        which="major",
        colors=cmap.colors[1],
        length=major_length,
        width=1,
        labelsize=fontsize,
    )
    ax2.tick_params(
        "y",
        which="minor",
        colors=cmap.colors[1],
        length=minor_length,
        width=1,
        labelsize=fontsize,
    )
    # Set y ticks to half multiples of pi
    ax2.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 1))

    def func(val, pos):
        if abs(1 - val / np.pi) < 1e-3:
            return f"$\pi$"
        elif val == 0:
            return "0"
        else:
            "{:.2g}$\pi$".format(val / np.pi)

    ax2.yaxis.set_major_formatter(FuncFormatter(func))

    # Adjust y-axis layout
    # Make y-axis of ax2 to the right of ax1
    ax2.spines["right"].set_position(("outward", 32))
    # Customise y-axis label position
    ax2.yaxis.set_label_coords(1.18, 0.5)

    ax3 = ax1.twinx()

    # Colour three lines
    ax3.plot(
        tidal_data.index[index_i:index_f],
        tidal_data["value_0"][index_i:index_f],
        color=cmap.colors[2],
        linestyle="-",
    )
    ax3.plot(
        tidal_data.index[index_i:index_f],
        tidal_data["value_1"][index_i:index_f],
        color=cmap.colors[2],
        linestyle="-",
    )
    ax3.plot(
        tidal_data.index[index_i:index_f],
        tidal_data["value_2"][index_i:index_f],
        color=cmap.colors[2],
        linestyle="-",
    )
    ax3.plot(
        tidal_data.index[index_i:index_f],
        tidal_data["value_3"][index_i:index_f],
        color=cmap.colors[2],
        linestyle="-",
    )
    # Increase ylim of ax3 by 0.1
    ax3.set_ylim(ax3.get_ylim()[0] * 1.1, ax3.get_ylim()[1] * 1.1)
    ax3.set_ylabel("$\\lambda_i$", color=cmap.colors[2], fontsize=fontsize)
    ax3.spines["right"].set_color(cmap.colors[2])
    ax3.spines["left"].set_color(cmap.colors[0])

    ax3.tick_params(
        "y",
        which="major",
        colors=cmap.colors[2],
        length=major_length,
        width=1,
        labelsize=fontsize,
    )
    ax3.tick_params(
        "y",
        which="minor",
        colors=cmap.colors[2],
        length=minor_length,
        width=1,
        labelsize=fontsize,
    )
    # ax3.ticklabel_format(axis="y", style="sci", scilimits=(0, 0)) # Scientific notation

    # Customise y-axis label position
    ax3.yaxis.set_label_coords(1.025, 0.5)

    # Draw reference line
    reference = params["tidal_reference"]
    if reference is not None:
        for ref in reference:
            ax3.axhline(y=ref, color="grey", linestyle="--")
            # Label numerical value at the end of the line
            ax3.text(
                x=ax3.get_xlim()[1] * 0.94,
                y=ref,
                s=f"{ref:.3f}",
                color="grey",
                ha="right",
                va="bottom",
                fontsize=fontsize - 2,
            )

    # Remove zero from ax3 y-axis
    old_ticks = ax3.get_yticks()
    # Get index of zero
    # zero_index = old_ticks.tolist().index(0)
    temp = old_ticks.tolist()
    # Find the index of the zero (or closest to zero)
    zero_index = min(range(len(temp)), key=lambda i: abs(temp[i] - 0))
    display(zero_index, temp[zero_index])
    # Set it to be invisible
    ax3.get_yticklabels()[zero_index].set_visible(False)

    # Add legend at top right
    if legend:
        # Add legend for whole plot
        fig.legend(
            ["Radial coordinate", "Polar angle", "Tidal eigenvalues"],
            loc="lower left",
            bbox_to_anchor=(0.06, 0.15),
        )

    # fig.tight_layout()

    if save:
        fig.savefig(
            give_filename(params, tidal=True).replace(".dat", ".png"),
            dpi=dpi,
        )

    if show:
        plt.show()

    plt.close()


def tidal_plot_alt(
    params,
    show=False,
    save=False,
    legend=True,
    dpi=1000,
    figsize=(10, 7),
    fontsize=18,
    major_length=3,
    minor_length=2,
    xaxis=None,
    plot_range=[0, -1],
    custom_cmap=None,
):
    filename = give_filename(params, tidal=False)
    data = load_geodesic(filename)

    filename = give_filename(params, tidal=True)
    tidal_data = load_tensor_cpp(filename)

    if not custom_cmap:
        custom_cmap = plt.get_cmap("Dark2")

    if plot_range == "auto":
        # Set start to when first r < 10*params["x"], searching from the beginning
        index_i = data[data["r"] < 6 * params["x"]].index[0]
        index_f = index_i + 2000
        print("Auto range:", index_i, index_f)
    else:
        index_i = plot_range[0]
        index_f = plot_range[1]

    if xaxis:
        xdata = data[xaxis][index_i:index_f]
    else:
        xdata = data["t"][index_i:index_f]

    # # Calculate how much xdata changes over the range and round to nearest power of 10
    # xdata_change = abs(xdata.iloc[-1] - xdata.iloc[0])
    # xdata_change = np.round(xdata_change, -int(np.floor(np.log10(xdata_change))))
    # # Change 200->100, 2000->1000, 20000->10000
    # xdata_change = int(xdata_change)
    # xdata_change = 10**int(np.log10(xdata_change))
    # print('xdata change:', xdata_change)

    # # Take away above this power
    # # e.g. 16234050 -> 4050 if power is 10^3
    # xdata = np.array([x % (xdata_change*100) for x in xdata])

    # Plot r, theta, and second column of tidal tensor on the same plot
    fig, ax1 = plt.subplots(figsize=figsize)

    # Set background colour to white
    fig.patch.set_facecolor("white")

    ax1.plot(
        xdata, data["r"][index_i:index_f], color=custom_cmap.colors[0], linestyle="--"
    )
    ax1.set_ylim(1, 1e2)
    ax1.set_yscale("log")  # Fix ax1 y-axis limits

    ax1.set_xlabel("$t$", fontsize=fontsize)
    ax1.set_ylabel("$r$", color=custom_cmap.colors[0], fontsize=fontsize)
    ax1.spines["left"].set_color(custom_cmap.colors[0])
    ax1.tick_params(
        "y",
        which="major",
        colors=custom_cmap.colors[0],
        length=major_length,
        width=1,
        labelsize=fontsize,
    )
    ax1.tick_params(
        "y",
        which="minor",
        colors=custom_cmap.colors[0],
        length=minor_length,
        width=1,
        labelsize=fontsize,
    )
    ax1.tick_params(
        "x", which="major", length=major_length, width=1, labelsize=fontsize
    )
    ax1.tick_params(
        "x", which="minor", length=minor_length, width=1, labelsize=fontsize
    )
    # Only tick even powers of 10
    ax1.yaxis.set_major_locator(LogLocator(base=100))
    ax1.minorticks_off()
    # Put x-label at the top
    # ax1.xaxis.set_label_position("top")
    # Customise y-axis label position
    ax1.yaxis.set_label_coords(-0.02, 0.5)
    
    # Hide offset text for x-axis
    plt.setp(ax1.get_xaxis().get_offset_text(), visible=False)
    
    # The below code hides the middle x tick and moves the axis label up slightly
    # # Customise x-axis label position
    # ax1.xaxis.set_label_coords(0.5, -0.02)
    # # print all x major ticks
    # print("x major ticks:", ax1.get_xticks())
    # # Set the middle one to be invisible
    # mid_index = len(ax1.get_xticks()) // 2
    # if mid_index < len(ax1.get_xticklabels()):
    #     ax1.get_xticklabels()[mid_index].set_visible(False)

    ax2 = ax1.twinx()

    ax2.plot(
        xdata,
        data["theta"][index_i:index_f],
        color=custom_cmap.colors[1],
        linestyle=":",
    )
    ax2.set_ylim(0, np.pi)
    ax2.set_ylabel("$\\theta$", color=custom_cmap.colors[1], fontsize=fontsize)
    ax2.spines["right"].set_color(custom_cmap.colors[1])
    ax2.spines["left"].set_color(custom_cmap.colors[0])
    ax2.tick_params(
        "y",
        which="major",
        colors=custom_cmap.colors[1],
        length=major_length,
        width=1,
        labelsize=fontsize,
    )
    ax2.tick_params(
        "y",
        which="minor",
        colors=custom_cmap.colors[1],
        length=minor_length,
        width=1,
        labelsize=fontsize,
    )
    # Set y ticks to half multiples of pi
    ax2.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 1))

    def func(val, pos):
        if abs(1 - val / np.pi) < 1e-3:
            return f"$\pi$"
        elif val == 0:
            return "0"
        else:
            "{:.2g}$\pi$".format(val / np.pi)

    ax2.yaxis.set_major_formatter(FuncFormatter(func))

    # Adjust y-axis layout
    # Make y-axis of ax2 to the right of ax1
    ax2.spines["right"].set_position(("outward", 32))
    # Customise y-axis label position
    ax2.yaxis.set_label_coords(1.18, 0.5)

    ax3 = ax1.twinx()

    # Colour three lines
    ax3.plot(
        xdata,
        tidal_data["value_0"][index_i:index_f],
        color=custom_cmap.colors[2],
        linestyle="-",
    )
    ax3.plot(
        xdata,
        tidal_data["value_1"][index_i:index_f],
        color=custom_cmap.colors[2],
        linestyle="-",
    )
    ax3.plot(
        xdata,
        tidal_data["value_2"][index_i:index_f],
        color=custom_cmap.colors[2],
        linestyle="-",
    )
    ax3.plot(
        xdata,
        tidal_data["value_3"][index_i:index_f],
        color=custom_cmap.colors[2],
        linestyle="-",
    )
    # Increase ylim of ax3 by 0.1
    ax3.set_ylim(ax3.get_ylim()[0] * 1.1, ax3.get_ylim()[1] * 1.1)
    ax3.set_ylabel("$\\lambda_i$", color=custom_cmap.colors[2], fontsize=fontsize)
    ax3.spines["right"].set_color(custom_cmap.colors[2])
    ax3.spines["left"].set_color(custom_cmap.colors[0])

    ax3.tick_params(
        "y",
        which="major",
        colors=custom_cmap.colors[2],
        length=major_length,
        width=1,
        labelsize=fontsize,
    )
    ax3.tick_params(
        "y",
        which="minor",
        colors=custom_cmap.colors[2],
        length=minor_length,
        width=1,
        labelsize=fontsize,
    )
    # ax3.ticklabel_format(axis="y", style="sci", scilimits=(0, 0)) # Scientific notation

    # Customise y-axis label position
    ax3.yaxis.set_label_coords(1.025, 0.5)

    # print(ax1.get_xlim())

    # Draw reference line
    reference = params["tidal_reference"]
    if reference is not None:
        for ref in reference:
            ax3.axhline(y=ref, color="grey", linestyle="--")
            # Label numerical value
            ax3.text(
                x=abs(ax1.get_xlim()[1] - ax1.get_xlim()[0]) * 0.45 + ax1.get_xlim()[0],
                y=ref,
                s=f"{ref:.3f}",
                color="grey",
                ha="right",
                va="bottom",
                fontsize=fontsize - 2,
            )

    # Remove zero from ax3 y-axis
    old_ticks = ax3.get_yticks()
    # Get index of zero
    # zero_index = old_ticks.tolist().index(0)
    temp = old_ticks.tolist()
    # Find the index of the zero (or closest to zero)
    zero_index = min(range(len(temp)), key=lambda i: abs(temp[i] - 0))
    display(zero_index, temp[zero_index])
    # Set it to be invisible
    ax3.get_yticklabels()[zero_index].set_visible(False)

    # Add legend at top right
    if legend:
        # Add legend for whole plot
        fig.legend(
            ["Radial coordinate", "Polar angle", "Tidal eigenvalues"],
            loc="lower left",
            bbox_to_anchor=(0.06, 0.15),
        )

    # fig.tight_layout()

    if save:
        fig.savefig(
            give_filename(params, tidal=True).replace(".dat", ".png"),
            dpi=dpi,
        )

    if show:
        plt.show()

    plt.close()
