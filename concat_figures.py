############################################################################################
# Copyright (C) Jędrzej Kubica, Nicolas Thierry-Mieg, 2026
#
# This file was written by Jędrzej Kubica and Nicolas Thierry-Mieg
# (CNRS, France) Nicolas.Thierry-Mieg@univ-grenoble-alpes.fr
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <https://www.gnu.org/licenses/>.
############################################################################################

import os
import sys
import math
import argparse
import pathlib
import logging

import matplotlib.pyplot
import matplotlib.image

# set up logger, using inherited config, in case we get called as a module
logger = logging.getLogger(__name__)

matplotlib.pyplot.set_loglevel('warning')
logging.getLogger("PIL").setLevel(logging.WARNING)


def main(figures, out_figure="fig_all.png"):
    logger.info(f"Concatenating {len(figures)} figures into one")
    labels = []
    if len(figures) > 25:
        raise Exception("Too many figures, cannot label with single letters (a), (b), ...")
    for i in range(len(figures)):
        labels.append(f"({chr(97 + i)})")  # (a), (b), (c), (d)...

    if len(figures) == 2:
        rows, cols = 1, 2
    else:
        rows = 2
        cols = max(1, math.ceil(len(figures) / rows))

    fig, axes = matplotlib.pyplot.subplots(rows, cols, figsize=(cols * 6, rows * 4.5))
    axes_flat = axes.ravel() if hasattr(axes, "ravel") else [axes]

    for ax, path, label in zip(axes_flat, figures, labels):
        img = matplotlib.image.imread(path)
        ax.imshow(img)
        ax.axis("off")

        ax.text(
            0.02, 0.95, label,
            transform=ax.transAxes,
            fontsize=10,
            fontweight="normal",
            va="top", ha="left"
        )

    for ax in axes_flat[len(figures):]:
        ax.axis("off")

    matplotlib.pyplot.subplots_adjust(wspace=0, hspace=0)
    matplotlib.pyplot.savefig(out_figure, bbox_inches="tight", dpi=1000)

    logger.info(f"Figure saved to {out_figure}")


if __name__ == "__main__":
    script_name = os.path.basename(sys.argv[0])
    # configure logging, sub-modules will inherit this config
    logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    # set up logger: we want script name rather than 'root'
    logger = logging.getLogger(script_name)

    parser = argparse.ArgumentParser(
        prog=script_name,
        description="""
        Concatenate multiple figures into one, with labels (a), (b), ... for each figure.
        The figures will be arranged in 2 rows, and as many columns as needed.
        """)
    parser.add_argument('--figures',
                        help="Paths (one after the other) to the figures to concatenate, in the order they should appear in the output figure",
                        nargs='+',
                        type=pathlib.Path,
                        required=True)
    parser.add_argument('--out',
                        help="Path to the output figure",
                        type=pathlib.Path,
                        default="fig_all.png")

    args = parser.parse_args()

    try:
        main(args.figures, out_figure=args.out)

    except Exception as e:
        # details on the issue should be in the exception name, print to stderr and die
        sys.stderr.write(f"ERROR in {script_name} : {repr(e)}\n")
        sys.exit(1)