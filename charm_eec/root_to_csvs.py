import ROOT
import numpy as np
import csv
import sys

ptbin = 5; # inclusive on last bin, inclusive on lower limit, exclusive on upper
pt_lo = [5, 10, 20, 40, 5]
pt_hi = [10, 20, 40, 60, 60]

etabin = 6; # inclusive on last bin, inclusive on lower limit, exclusive on upper
eta_lo = [-3.5, -1, 1, -3.5, -1, 0]
eta_hi = [-1, 1, 3.5, 3.5, 0, 1]

# assumes all histograms in hists have the same number of bins, bin centers, and bin widths
def hists_to_csv(outfile_name, hists):
  # Open output CSV file for writing
  with open(outfile_name, mode="w") as output_file:

    # Create CSV writer object
    writer = csv.writer(output_file)

    nbins = hists[0].GetNbinsX()

    # Loop over bins and write row to CSV file
    for i in range(1, nbins+1):
      bin_center = hists[0].GetBinCenter(i)
      bin_width = hists[0].GetBinWidth(i)

      this_row = [bin_center, bin_width]
      for ihist in range(len(hists)):
        this_row.append(hists[i].GetBinContent(i))
        this_row.append(hists[i].GetBinError(i))

      writer.writerow(this_row)
    

def root_to_csvs(infile="merged.root", outdir="./"):
  infile = sys.argv[1]
  outdir = sys.argv[2]

  fin = ROOT.TFile(infile)
  
  # convert each histogram to csv

  # create EEC overlay histograms: overlays EEC across ptbin, on file per eta bin
  for ieta in range(3):
    hists_this_eta = []
    for ipt in range(3):
      hist = fin.Get("h1d_jet_eec", ieta, "_", ipt)
      hists_this_eta.append(hist)
    output_path = outdir + "eec_overlay_" + str(ieta) + ".csv"
    hists_to_csv(output_path, hists_this_eta)

if __name__ == "__main__":
  root_to_csvs()