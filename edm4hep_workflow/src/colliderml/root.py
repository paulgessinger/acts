import numpy


class TH1:
    def __init__(self, th1_tefficiency, xrange=None):
        import ROOT

        if isinstance(th1_tefficiency, ROOT.TEfficiency):
            th1 = th1_tefficiency.GetTotalHistogram()
        else:
            th1 = th1_tefficiency

        bins = list(range(1, th1.GetNbinsX() + 1))

        if xrange is not None:
            bins = [
                i
                for i in bins
                if th1.GetBinCenter(i) >= xrange[0] and th1.GetBinCenter(i) <= xrange[1]
            ]

        self.x = numpy.array([th1.GetBinCenter(i) for i in bins])

        self.x_lo = numpy.array([th1.GetBinLowEdge(i) for i in bins])
        self.x_width = numpy.array([th1.GetBinWidth(i) for i in bins])
        self.x_hi = numpy.add(self.x_lo, self.x_width)
        self.x_err_lo = numpy.subtract(self.x, self.x_lo)
        self.x_err_hi = numpy.subtract(self.x_hi, self.x)

        try:
            self.y = numpy.array([th1_tefficiency.GetEfficiency(i) for i in bins])
            self.y_err_lo = numpy.array(
                [th1_tefficiency.GetEfficiencyErrorLow(i) for i in bins]
            )
            self.y_err_hi = numpy.array(
                [th1_tefficiency.GetEfficiencyErrorUp(i) for i in bins]
            )
        except Exception:
            self.y = numpy.array([th1_tefficiency.GetBinContent(i) for i in bins])
            self.y_err_lo = numpy.array([th1_tefficiency.GetBinError(i) for i in bins])
            self.y_err_hi = numpy.array([th1_tefficiency.GetBinError(i) for i in bins])

    def errorbar(self, ax, **errorbar_kwargs):
        ax.errorbar(
            self.x,
            self.y,
            yerr=(self.y_err_lo, self.y_err_hi),
            xerr=(self.x_err_lo, self.x_err_hi),
            **errorbar_kwargs,
        )
        return ax

    def step(self, ax, **step_kwargs):
        ax.step(self.x_hi, self.y, **step_kwargs)
        return ax

    def bar(self, ax, **bar_kwargs):
        ax.bar(self.x, height=self.y, yerr=(self.y_err_lo, self.y_err_hi), **bar_kwargs)
        return ax
