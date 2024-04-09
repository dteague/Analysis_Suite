#!/usr/bin/env python3
from scipy.interpolate import griddata
import numpy as np
import argparse
import uproot
import boost_histogram as bh

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

import analysis_suite.commons.user as user
from analysis_suite.combine.combine_wrapper import runCombine

class Data2D:
    def __init__(self, data_file, corr_file, n_points=100, x_range=[0.,3.], y_range=[0.,3.]):
        with uproot.open(data_file) as f:
            arr = f["limit"].arrays()
            self.y, self.x = arr.r_sig, arr.r_other
            self.dnll = arr.deltaNLL
        self.corr = np.nan
        try:
            if "Hesse" in corr_file.name:
                with uproot.open(corr_file) as f:
                    self.corr = f['h_correlation'].to_boost()[bh.loc("r_sig"), bh.loc("r_other")]
            else:
                f =  ROOT.TFile.Open(str(corr_file))
                fit = f.Get("fit_mdf")
                self.corr = fit.correlation("r_sig", "r_other")
                f.Close()
        except:
            pass
        self.n_points = n_points
        self.xlo, self.xhi = x_range
        self.ylo, self.yhi = y_range

    def get_grid(self):
        points = np.array([self.x, self.y]).transpose()
        grid_x, grid_y = np.mgrid[self.xlo:self.xhi:self.n_points*1j, self.ylo:self.yhi:self.n_points*1j]
        grid_vals = griddata(points, self.dnll, (grid_x,grid_y), "cubic")

        # Remove NANS
        grid_x = grid_x[grid_vals==grid_vals]
        grid_y = grid_y[grid_vals==grid_vals]
        grid_vals = grid_vals[grid_vals==grid_vals]
        grid_vals = np.where(grid_vals < 0, 0, grid_vals)
        return (grid_x, grid_y, grid_vals)


def run_combine(workdir, cardname, year, other, sig='ttt', sig_up=3, sig_down=0, other_up=3, other_down=0, hesse=False, nosyst=False, debug=False, skip=False):
    corr_command = '--robustHesse 1' if hesse else "--robustFit 1"
    # extra_command = "--freezeParameters allConstrainedNuisances"
    corr_filename = "robustHesse" if hesse else "multidimfit"
    out_corr, out_plot = "corr", "plot"

    if nosyst:
        out_corr += "_nosyst"
        out_plot += "_nosyst"

    if not skip:
        runCombine.work_dir = workdir
        runCombine(f'text2workspace.py {cardname}.txt -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
                     --PO "map=.*/{sig}:r_sig[1,0,200]" --PO "map=.*/{other}:r_other[1,0,200]" -o {cardname}.root', output=debug)
        runCombine(f'combine -M MultiDimFit {cardname}.root -m 125 -n .scan2d.{out_plot}.{year}  --points 800 --algo grid \
                     --cminDefaultMinimizerStrategy 0 -t -1 -P r_sig -P r_other --setParameterRanges r_sig={sig_down},{sig_up}:r_other={other_down},{other_up} \
                     --expectSignal=1', output=debug)
        runCombine(f'combine -M MultiDimFit {cardname}.root -m 125 -n .scan2d.{out_corr}.{year} --expectSignal=1 \
                     --cminDefaultMinimizerStrategy 0 -P r_sig -P r_other --setParameterRanges r_sig={sig_down},{sig_up}:r_other={other_down},{other_up}\
                     -t -1 {corr_command} --saveFitResult\
                     --stepSize=0.01 --setRobustFitTolerance=1000 --setCrossingTolerance=0.001 ', output=debug)
        #### single
        # runCombine(f'text2workspace.py {cardname}_single.txt -m 125 -o {cardname}_single.root', output=debug)
        # runCombine(f'combine -M Significance {cardname}_single.root -t -1 --expectSignal 1 --toysFrequentist ', output=debug)
    return workdir/f'higgsCombine.scan2d.{out_plot}.{year}.MultiDimFit.mH125.root', workdir/f'{corr_filename}.scan2d.{out_corr}.{year}.root'

# Make plot with Root
def make_2d_plot(outfile, data, xaxis_name="", yaxis_name=""):
    grid_x, grid_y, grid_vals = data.get_grid()
    n_bins = 40

    # Define Profile2D histogram
    h2D = ROOT.TProfile2D("h","h",n_bins, data.xlo, data.xhi, n_bins, data.ylo, data.yhi)
    for i in range(len(grid_vals)):
        # Factor of 2 comes from 2*NLL
        h2D.Fill( grid_x[i], grid_y[i], 2*grid_vals[i] )
        # h2D.Fill( grid_x[i], grid_y[i], grid_vals[i] )

    # Set up canvas
    canv = ROOT.TCanvas("canv","canv",600,600)
    canv.SetTickx()
    canv.SetTicky()
    canv.SetLeftMargin(0.115)
    canv.SetBottomMargin(0.115)

    xhigh = data.xhi - (data.xhi-data.xlo)/n_bins
    yhigh = data.yhi - (data.yhi-data.ylo)/n_bins

    # Set histogram properties
    h2D.SetContour(999)
    h2D.SetTitle("")
    h2D.GetXaxis().SetTitle(xaxis_name)
    h2D.GetXaxis().SetTitleSize(0.05)
    h2D.GetXaxis().SetTitleOffset(0.9)
    h2D.GetXaxis().SetRangeUser(data.xlo, xhigh)

    h2D.GetYaxis().SetTitle(yaxis_name)
    h2D.GetYaxis().SetTitleSize(0.05)
    h2D.GetYaxis().SetTitleOffset(0.9)
    h2D.GetYaxis().SetRangeUser(data.ylo, yhigh)

    h2D.GetZaxis().SetTitle("-2 #Delta ln L")
    h2D.GetZaxis().SetTitleSize(0.05)
    h2D.GetZaxis().SetTitleOffset(0.8)
    h2D.GetZaxis().SetRangeUser(-0.01, 25)

    # h2D.SetMaximum(25)

    # Make confidence interval contours
    c68, c95 = h2D.Clone(), h2D.Clone()
    c68.SetContour(2)
    c68.SetContourLevel(1,2.3)
    c68.SetLineWidth(3)
    c68.SetLineColor(ROOT.kBlack)
    c95.SetContour(2)
    c95.SetContourLevel(1,5.99)
    c95.SetLineWidth(3)
    c95.SetLineStyle(2)
    c95.SetLineColor(ROOT.kBlack)

    # Draw histogram and contours
    h2D.Draw("COLZ")

    # Draw lines for SM point
    vline = ROOT.TLine(1, data.ylo, 1, yhigh)
    vline.SetLineColorAlpha(ROOT.kGray,0.5)
    vline.Draw("Same")
    hline = ROOT.TLine(data.xlo, 1, xhigh, 1)
    hline.SetLineColorAlpha(ROOT.kGray,0.5)
    hline.Draw("Same")

    # Draw contours
    c68.Draw("cont3same")
    c95.Draw("cont3same")

    # Make best fit and sm points
    gBF = ROOT.TGraph()
    gBF.SetPoint(0,grid_x[np.argmin(grid_vals)],grid_y[np.argmin(grid_vals)])
    gBF.SetMarkerStyle(34)
    gBF.SetMarkerSize(2)
    gBF.SetMarkerColor(ROOT.kBlack)
    gBF.Draw("P")

    gSM = ROOT.TGraph()
    gSM.SetPoint(0,1,1)
    gSM.SetMarkerStyle(33)
    gSM.SetMarkerSize(2)
    gSM.SetMarkerColor(ROOT.kRed)
    gSM.Draw("P")

    # Add legend
    leg = ROOT.TLegend(0.6,0.67,0.8,0.87)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.AddEntry(gBF,  "Best fit", "P" )
    leg.AddEntry(c68, "1#sigma CL" , "L" )
    leg.AddEntry(c95, "2#sigma CL" , "L" )
    leg.AddEntry(gSM,  "SM"     , "P" )
    leg.Draw()

    corr = str(round(data.corr, 3))
    corr_text = ROOT.TLatex()
    corr_text.DrawLatex(0.05*(xhigh-data.xlo)+data.xlo, 0.9*(yhigh-data.ylo)+data.ylo, r'#rho = '+corr)
    corr_text.Draw()
    print(corr)
    canv.Update()
    canv.SaveAs(str(outfile))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="main", description="")
    parser.add_argument("-d", "--workdir", required=True, type=lambda x : user.workspace_area/x/"combine",
                        help="Working Directory")
    parser.add_argument('-t', '--extra_text')
    parser.add_argument("-y", "--years", required=True,
                        type=lambda x : [i.strip() for i in x.split(',')], help="Year to use")
    parser.add_argument('-o', '--other', default='4top')
    parser.add_argument('-s', '--signal', default='ttt')
    parser.add_argument('-ns', '--no_syst', action="store_true")
    parser.add_argument("--skip", action="store_true")
    parser.add_argument("--debug", action='store_true')
    args = parser.parse_args()

    extra = args.other if args.extra_text is None else args.extra_text
    sig_up = 50
    other_up = 3
    other_down = 0
    sig_down = 0
    syst_name = "_nosyst" if args.no_syst else ""

    # workdir = args.workdir/extra
    workdir = user.analysis_area/f'corr_test/final_Dec_{args.extra_text}/combine'

    for year in args.years:
        cardname = f'final_{year}{syst_name}_card'
        # cardname = f'signal_{year}_Dilepton_card'
        data_file, corr_file = run_combine(workdir, cardname, year, args.other, sig=args.signal, sig_up=sig_up, other_down=other_down, sig_down=sig_down,
                                           other_up=other_up, debug=args.debug, hesse=False, nosyst=args.no_syst, skip=args.skip)
        data = Data2D(data_file, corr_file, y_range=[sig_down, sig_up], x_range=[other_down, other_up])
        make_2d_plot(workdir/f"correlation_{args.signal}{syst_name}_{year}.png", data, xaxis_name=f"r_{args.other}", yaxis_name=f'r_{args.signal}')


def run_single(workdir, cardname, year, signal, other, no_syst=False, sig_up=50, sig_down=0, other_up=3, other_down=0, extra=""):
    syst_name = "_nosyst" if no_syst else ""
    debug=False
    if '.txt' in cardname:
        cardname = cardname[:cardname.index('.txt')]
    data_file, corr_file = run_combine(workdir, cardname, year, other, sig=signal, sig_up=sig_up, other_down=other_down, sig_down=sig_down,
                                       other_up=other_up, hesse=False, nosyst=no_syst, debug=debug)
    data = Data2D(data_file, corr_file, y_range=[sig_down, sig_up], x_range=[other_down, other_up])
    make_2d_plot(workdir/f"correlation_{signal}{syst_name}{extra}_{year}.png", data, xaxis_name=f"r_{other}", yaxis_name=f'r_{signal}')
