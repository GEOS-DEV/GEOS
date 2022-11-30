
import os
import warnings
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors


class HighResPlot():
    def __init__(self):
        self.handle = plt.figure()
        self.isAdjusted = False
        self.colorspace = 'RGB'
        self.compress = False
        self.output_format = 'png'
        self.font_weight = 'normal'
        self.axis_font_size = 8
        self.legend_font_size = 8
        self.size = (3.5, 2.625)
        self.dpi = 200
        self.apply_tight_layout = True
        self.disable_antialiasing = False
        self.use_imagemagick = False
        self.apply_font_size()

    def __del__(self):
        try:
            self.handle.clf()
            plt.close(self.handle.number)
        except:
            pass

    def set_active(self):
        plt.figure(self.handle.number)

    def reset(self):
        tmp = self.handle.number
        self.set_active()
        self.handle.clf()
        del self.handle
        self.handle = plt.figure(tmp)
        self.apply_font_size()
        self.isAdjusted = False

    def apply_font_size(self):
        self.font = {'weight': self.font_weight,
                     'size': self.axis_font_size}
        plt.rc('font', **self.font)

    def apply_format(self):
        self.set_active()
        self.handle.set_size_inches(self.size)
        self.apply_font_size()

        for ax in self.handle.get_axes():
            tmp = ax.get_legend()
            if tmp:
                ax.legend(loc=tmp._loc,
                          fontsize=self.legend_font_size)

        if self.apply_tight_layout:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                self.handle.tight_layout(pad=1.5)
                self.apply_tight_layout = False

    def save(self, fname):
        self.apply_format()
        if plt.isinteractive():
            plt.draw()
            plt.pause(0.1)

        if self.use_imagemagick:
            # Matplotlib's non-pdf layout is sometimes odd,
            # so use the imagemagick convert as a back-up option
            plt.savefig('%s.pdf' % (fname),
                        format='pdf',
                        dpi=self.dpi)

            if (self.output_format != 'pdf'):
                cmd = 'convert -density %i %s.pdf' % (self.dpi, fname)

                if (self.colorspace == 'CMYK'):
                    cmd += ' -colorspace CMYK'
                if (self.colorspace == 'gray'):
                    cmd += ' -colorspace gray'
                if (self.compress):
                    cmd += ' -compress lzw'
                if (self.disable_antialiasing):
                    cmd += ' +antialias'

                cmd += ' %s.%s' % (fname, self.output_format)
                os.system(cmd)
                os.system('rm %s.pdf' % fname)
        else:
            plt.savefig('%s.%s' % (fname, self.output_format),
                        format=self.output_format,
                        dpi=self.dpi)


def get_modified_jet_colormap():
    # Add an updated colormap
    cdict = cm.get_cmap('jet').__dict__['_segmentdata'].copy()
    for k in cdict.keys():
        tmp_seq = list(cdict[k])
        tmp_final = list(tmp_seq[0])
        tmp_final[1] = 1.0
        tmp_final[2] = 1.0
        tmp_seq[0] = tuple(tmp_final)
        cdict[k] = tuple(tmp_seq)
    return mcolors.LinearSegmentedColormap('jet_mod', cdict)


def get_periodic_line_style(ii):
    col = ['k', 'b', 'c', 'g', 'y', 'r', 'm']
    mark = ['-', '--', ':', '-.', 'o', '*', '^', '+']

    jj = ii % (len(col) * len(mark))
    style = col[jj % len(col)] + mark[int(jj/len(col))]

    return style

