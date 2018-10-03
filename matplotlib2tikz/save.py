# -*- coding: utf-8 -*-
#
from __future__ import print_function

import codecs
import os
import matplotlib as mpl
import six

# For DateFile
from itertools import zip_longest
from collections import OrderedDict
import numpy as np

from . import axes
from . import legend
from . import line2d
from . import image as img
from . import quadmesh as qmsh
from . import path
from . import patch
from . import text as mytext
from .__about__ import __version__


def get_tikz_code(
        filepath,
        figure='gcf',
        figurewidth=None,
        figureheight=None,
        textsize=10.0,
        tex_relative_path_to_data=None,
        strict=False,
        wrap=True,
        extra_axis_parameters=None,
        extra_tikzpicture_parameters=None,
        dpi=None,
        show_info=True,
        manual_legend=False,
        figlabel='',
        data_file = None
        ):
    '''Main function. Here, the recursion into the image starts and the
    contents are picked up. The actual file gets written in this routine.

    :param figure: either a Figure object or 'gcf' (default).

    :param filepath: The file to which the TikZ output will be written.
    :type filepath: str

    :param figurewidth: If not ``None``, this will be used as figure width
                        within the TikZ/PGFPlots output. If ``figureheight``
                        is not given, ``matplotlib2tikz`` will try to preserve
                        the original width/height ratio.
                        Note that ``figurewidth`` can be a string literal,
                        such as ``'\\figurewidth'``.
    :type figurewidth: str

    :param figureheight: If not ``None``, this will be used as figure height
                         within the TikZ/PGFPlots output. If ``figurewidth`` is
                         not given, ``matplotlib2tikz`` will try to preserve
                         the original width/height ratio.  Note that
                         ``figurewidth`` can be a string literal, such as
                         ``'\\figureheight'``.
    :type figureheight: str

    :param textsize: The text size (in pt) that the target latex document is
                     using.  Default is 10.0.
    :type textsize: float

    :param tex_relative_path_to_data: In some cases, the TikZ file will have to
                                      refer to another file, e.g., a PNG for
                                      image plots. When ``\\input`` into a
                                      regular LaTeX document, the additional
                                      file is looked for in a folder relative
                                      to the LaTeX file, not the TikZ file.
                                      This arguments optionally sets the
                                      relative path from the LaTeX file to the
                                      data.
    :type tex_relative_path_to_data: str

    :param strict: Whether or not to strictly stick to matplotlib's appearance.
                   This influences, for example, whether tick marks are set
                   exactly as in the matplotlib plot, or if TikZ/PGFPlots
                   can decide where to put the ticks.
    :type strict: bool

    :param wrap: Whether ``'\\begin{tikzpicture}'`` and
                 ``'\\end{tikzpicture}'`` will be written. One might need to
                 provide custom arguments to the environment (eg. scale= etc.).
                 Default is ``True``.
    :type wrap: bool

    :param extra_axis_parameters: Extra axis options to be passed (as a set)
                                    to pgfplots. Default is ``None``.
    :type extra_axis_parameters: a set of strings for the pfgplots axes.

    :param extra_tikzpicture_parameters: Extra tikzpicture options to be passed
                                        (as a set) to pgfplots.

    :type extra_tikzpicture_parameters: a set of strings for the pfgplots
                                        tikzpicture.

    :param dpi: The resolution in dots per inch of the rendered image in case
                of QuadMesh plots. If ``None`` it will default to the value
                ``savefig.dpi`` from matplotlib.rcParams. Default is ``None``.
    :type dpi: int

    :returns: None


    The following optional attributes of matplotlib's objects are recognized
    and handled:

     - axes.Axes._matplotlib2tikz_anchors
       This attribute can be set to a list of ((x,y), anchor_name) tuples.
       Invisible nodes at the respective location will be created which  can be
       referenced from outside the axis environment.
    '''
    # not as default value because gcf() would be evaluated at import time
    if figure == 'gcf':
        figure = mpl.pyplot.gcf()
    data = {}
    data['manual-legend'] = manual_legend
    data['figlabel'] = figlabel
    data['fwidth'] = figurewidth
    data['fheight'] = figureheight
    data['rel data path'] = tex_relative_path_to_data
    data['output dir'] = os.path.dirname(filepath)
    data['base name'] = os.path.splitext(os.path.basename(filepath))[0]
    data['strict'] = strict
    data['tikz libs'] = set()
    data['pgfplots libs'] = set()
    data['font size'] = textsize
    data['custom colors'] = {}
    data['extra tikzpicture parameters'] = extra_tikzpicture_parameters
    # rectangle_legends is used to keep track of which rectangles have already
    # had \addlegendimage added. There should be only one \addlegenimage per
    # bar chart data series.
    data['rectangle_legends'] = set()
    if extra_axis_parameters:
        data['extra axis options [base]'] = extra_axis_parameters.copy()
    else:
        data['extra axis options [base]'] = set()

    if dpi:
        data['dpi'] = dpi
    else:
        savefig_dpi = mpl.rcParams['savefig.dpi']
        data['dpi'] = \
            savefig_dpi if isinstance(savefig_dpi, int) \
            else mpl.rcParams['figure.dpi']

    # gather the file content
    data, content = _recurse(data, figure, [], data_file)

    # Check if there is still an open groupplot environment. This occurs if not
    # all of the group plot slots are used.
    if 'is_in_groupplot_env' in data and data['is_in_groupplot_env']:
        content.extend('\\end{groupplot}\n\n')

    disclaimer = 'This file was created by matplotlib2tikz v%s.' % __version__

    # write disclaimer to the file header
    code = ''''''
    code += _tex_comment(disclaimer)

    # write the contents
    if wrap:
        code += '\\begin{tikzpicture}\n\n'
        if extra_tikzpicture_parameters:
            code += ',\n'.join(data['extra tikzpicture parameters'])
            code += '\n'

    coldefs = _get_color_definitions(data)
    if coldefs:
        code += '\n'.join(coldefs)
        code += '\n\n'

    code += ''.join(content)

    if wrap:
        code += '\\end{tikzpicture}'

    # print message about necessary pgfplot libs to command line
    if show_info:
        _print_pgfplot_libs_message(data)
    return code


def save(*args, **kwargs):
    '''Same as `get_tikz_code()`, but actually saves the code to a file.
    '''
    code = get_tikz_code(*args, **kwargs)

    encoding = kwargs['encoding'] if 'encoding' in kwargs else None
    file_handle = codecs.open(args[0], 'w', encoding)
    try:
        file_handle.write(code)
    except UnicodeEncodeError:
        # We're probably using Python 2, so treat unicode explicitly
        file_handle.write(six.text_type(code).encode('utf-8'))
    file_handle.close()
    return


def _tex_comment(comment):
    '''Prepends each line in string with the LaTeX comment key, '%'.
    '''
    return '% ' + str.replace(comment, '\n', '\n% ') + '\n'


def _get_color_definitions(data):
    '''Returns the list of custom color definitions for the TikZ file.
    '''
    definitions = []
    for name, rgb in data['custom colors'].items():
        definitions.append('\\definecolor{%s}{rgb}{%.15g,%.15g,%.15g}'
                           % (name, rgb[0], rgb[1], rgb[2])
                           )
    return definitions


def _print_pgfplot_libs_message(data):
    '''Prints message to screen indicating the use of PGFPlots and its
    libraries.'''
    pgfplotslibs = ','.join(list(data['pgfplots libs']))
    tikzlibs = ','.join(list(data['tikz libs']))

    print('=========================================================')
    print('Please add the following lines to your LaTeX preamble:\n')
    print('\\usepackage[utf8]{inputenc}')
    print('\\usepackage{fontspec}'
          ' % This line only for XeLaTeX and LuaLaTeX')
    print('\\usepackage{pgfplots}')
    if tikzlibs:
        print('\\usetikzlibrary{' + tikzlibs + '}')
    if pgfplotslibs:
        print('\\usepgfplotslibrary{' + pgfplotslibs + '}')
    print('=========================================================')
    return


class _ContentManager(object):
    """ Basic Content Manager for matplotlib2tikz

    This manager uses a dictionary to map z-order to an array of content
    to be drawn at the z-order.
    """
    def __init__(self):
        self._content = dict()

    def extend(self, content, zorder):
        """ Extends with a list and a z-order
        """
        if zorder not in self._content:
            self._content[zorder] = []
        self._content[zorder].extend(content)

    def flatten(self):
        content_out = []
        all_z = sorted(self._content.keys())
        for z in all_z:
            content_out.extend(self._content[z])
        return content_out


class DataFile(object):
    def __init__(self, filename):
        self.columns = OrderedDict()
        self.filename = filename

    def append(self, col_type, column, rel_tol=1e-09):
        cmp_eq = lambda a, b: np.abs(a-b) <= rel_tol * np.maximum(np.abs(a), np.abs(b))

        i = 0
        ac = np.array(column)
        for k, v in self.columns.items():
            if len(k) == len(col_type) or k[len(col_type):].isdigit():
                i += 1
            if len(v) == len(column) and np.all(cmp_eq(v, ac)):
                return k   # This column has already been added

        # New column
        key = col_type + (str(i) if i > 0 else '')
        self.columns[key] = ac
        return key

    def write(self, filename, transpose=False):
        keys = self.columns.keys()
        vals = self.columns.values()
        # Longest vector first
        # ii = np.argsort([len(v) for v in vals])[::-1]
        # keys = [keys[i] for i in ii]
        # vals = [vals[i] for i in ii]
        if not transpose:
            with open(filename, 'w') as f:
                f.write(" ".join(keys) + "\n")
                for row in zip_longest(*vals, fillvalue=np.nan):
                    f.write(" ".join(["%.15g" % v for v in row if v is not None]) + "\n")
        else:
            max_count = np.max([len(v_row) for v_row in vals])
            with open(filename, 'w') as f:
                for k_row, v_row in zip(keys, vals):
                    f.write(k_row + " ")
                    f.write(" ".join(["%.15g" % v for v in v_row if v is not
                                      None]))
                    if len(v_row) < max_count:
                        f.write(" " + " ".join(["nan"] * (max_count - len(v_row))))
                    f.write("\n")


def _recurse(data, obj, saved_objs, data_file):
    '''Iterates over all children of the current object, gathers the contents
    contributing to the resulting PGFPlots file, and returns those.
    '''
    content = _ContentManager()
    for child in obj.get_children():
        if child in saved_objs:
            continue
        saved_objs.append(child)

        if isinstance(child, mpl.axes.Axes):
            # Reset 'extra axis parameters' for every new Axes environment.
            data['extra axis options'] = \
                data['extra axis options [base]'].copy()
            ax = axes.Axes(data, child)
            if not ax.is_colorbar:
                # Run through the child objects, gather the content.

                # We first need to handle errorbars since they
                # are saved in child.containers
                container_content = _ContentManager()
                for con in child.containers:
                    if type(con) == mpl.container.ErrorbarContainer:
                        data, cont = line2d.draw_errorbar2d(data, con, data_file)
                        container_content.extend(cont, child.get_zorder())
                        saved_objs.extend(con.get_children())

                # Then recurse on the axes
                data, children_content = _recurse(data, child,
                                                  saved_objs, data_file)
                # add extra axis options from children
                if data['extra axis options']:
                    ax.axis_options.extend(data['extra axis options'])

                manual_legend = 'manual-legend' in data and data['manual-legend']

                leg, cont_leg = None, None
                if manual_legend:
                    leg = child.get_legend()
                    if leg is not None:
                        data, cont_leg = legend.draw_legend(data, leg)

                # populate content
                content.extend(
                    ax.get_begin_code() +
                    container_content.flatten() +
                    children_content +
                    (data['coordinates'] if 'coordinates' in data else []) +
                    [ax.get_end_code(data)], 0)

                if cont_leg is not None:
                    content.extend(cont_leg, leg.get_zorder())

        elif isinstance(child, mpl.lines.Line2D):
            data, cont = line2d.draw_line2d(data, child, data_file)
            content.extend(cont, child.get_zorder())
        elif isinstance(child, mpl.image.AxesImage):
            data, cont = img.draw_image(data, child)
            content.extend(cont, child.get_zorder())
            # # Really necessary?
            # data, children_content = _recurse(data, child, saved_objs, data_file)
            # content.extend(children_content)
        elif isinstance(child, mpl.patches.Patch):
            data, cont = patch.draw_patch(data, child)
            content.extend(cont, child.get_zorder())
        elif isinstance(
                child,
                (
                    mpl.collections.PatchCollection,
                    mpl.collections.PolyCollection
                )
                ):
            data, cont = patch.draw_patchcollection(data, child)
            content.extend(cont, child.get_zorder())
        elif isinstance(child, mpl.collections.PathCollection):
            data, cont = path.draw_pathcollection(data, child, data_file)
            content.extend(cont, child.get_zorder())
        elif isinstance(child, mpl.collections.LineCollection):
            data, cont = line2d.draw_linecollection(data, child)
            content.extend(cont, child.get_zorder())
        elif isinstance(child, mpl.collections.QuadMesh):
            data, cont = qmsh.draw_quadmesh(data, child)
            content.extend(cont, child.get_zorder())
        elif isinstance(child, mpl.legend.Legend):
            if 'manual-legend' in data and data['manual-legend']:
                pass
            else:
                data, cont = legend.draw_legend(data, child)
                content.extend(cont, child.get_zorder())
        elif isinstance(
                child,
                (
                    mpl.axis.XAxis, mpl.axis.YAxis,
                    mpl.spines.Spine, mpl.text.Text
                )):
            pass
        else:
            print('matplotlib2tikz: Don''t know how to handle object ''%s''.' %
                  type(child))
    # XXX: This is ugly
    if isinstance(obj, (mpl.axes.Subplot, mpl.figure.Figure)):
        for text in obj.texts:
            data, cont = mytext.draw_text(data, text)
            content.extend(cont, text.get_zorder())
    return data, content.flatten()
