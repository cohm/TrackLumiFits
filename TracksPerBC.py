# This file contains code to construct a function describing the distribution of the expected number of measurements per bunch crossing in a luminosity measurement, given a known distribution of measurements per interaction

import ROOT, sys

wait = False

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)

# returns a histogram representing probabilities for discrete values specified by the provided function
def getPdfHistogramFromFunction(fPdf, draw = False, canvas = None):
    if draw:
        fPdf.Draw()
        if wait:
            ROOT.gPad.WaitPrimitive() # stop here to show the plot
    # the function is assumed to be defined on a correct range for integers, e.g. [-0.5, 19.5] to give probabilities for integer values in the center of bins
    range = fPdf.GetMinimumX() - fPdf.GetMaximumX()
    print("Function is defined on [%f, %f] with a range of %f" % (fPdf.GetMinimumX(), fPdf.GetMaximumX(), range))
    if not range.is_integer():
        print("ERROR: function not defined correctly, needs to be evaluated for integer values")
        sys.exit(1)
    f.SetNpx(int(range)) # discretize
    if draw:
        fPdf.Draw("HISTSAME")
        if wait:
            ROOT.gPad.WaitPrimitive() # stop here to show the plot
    pdfHisto = fPdf.GetHistogram()
    return pdfHisto

# returns the PDF of a sum of two random variables described by the two provided PDFs in the histograms
def getPdfOfSum(h1, h2, n):
    # the sum is the convolution of them, bin 1 has entry for 0 tracks, etc
    hSum = h1.Clone() # create new with same range etc
    hSum.Reset()
    hSum.SetName("Measurements for %d interactions" % n)
    hSum.SetTitle(hSum.GetName())
    print("Will now calculate the track-multiplicity probabilities for %d interactions" % n)
    # loop over the track multiplicities (bin 1 is for 0 tracks)
    for nTrk in range(1, hSum.GetNbinsX()+1):
        #print("  nTrk: %d" % nTrk)
        # for each track multiplicity, add upp all combinations thet contribute
        nTrkProb = 0
        for nTrk1 in range(0, nTrk+1):
            nTrk2 = nTrk-nTrk1
            #print("    h1(%d) && h2(%d)" % (nTrk1, nTrk2))
            p1 = h1.GetBinContent(nTrk1+1)
            p2 = h2.GetBinContent(nTrk2+1)
            #print("      --> %f * %f = %f" % (p1, p2, p1*p2))
            nTrkProb += p1*p2
            #if p1 == 0:
            #    break # no need to keep looping
        #print("  --> nTrkProb(%d) = %f" % (nTrk, nTrkProb))
        hSum.SetBinContent(nTrk+1, nTrkProb)
    return hSum

# generate histograms for the PDFs of the number of expected measurements for a (range of) fixed number of interactions
def makeNInteractionHistos(perIntHisto, maxInteractions, maxMeasurements, draw = False, canvas = None):
    print("Will now make pdfs for number of measurements for a fixed number of interactions, from 0 to %d" % maxInteractions)
    # dict to populate, one entry per fixed number of interactions
    pdfs = {}
    hNameTemplate = "Track multiplicity probabilities for fixed nInt"
    # the first two we need to fill manually, then we can construct the rest with a loop
    # zero interactions is trivial
    pdfs[0] = ROOT.TH1D(hNameTemplate+(" for %d interactions" % 0), hNameTemplate, maxMeasurements+1, -0.5, maxMeasurements+0.5)
    # one interaction is the per-interaction pdf, but we need to make sure we use the requested track-multiplicity range
    pdfs[1] = ROOT.TH1D(hNameTemplate+(" for %d interactions" % 1), hNameTemplate, maxMeasurements+1, -0.5, maxMeasurements+0.5)
    # fill the bins with the contents from the per-interaction histo
    for b in range(1, perIntHisto.GetNbinsX()+1):
        pdfs[1].SetBinContent(b, perIntHisto.GetBinContent(b))
    pdfs[1].Scale(1/pdfs[1].Integral()) # a pdf is normalized
    pdfs[1].SetLineColor(1)
    if draw:
        pdfs[0].Draw()
        if wait:
            ROOT.gPad.WaitPrimitive()
        pdfs[1].Draw("HISTSAME")
        if wait:
            ROOT.gPad.WaitPrimitive()
        ROOT.gPad.SetLogy(1)
        if wait:
            ROOT.gPad.WaitPrimitive()

    # now create the pdfs for fixed number of interactions n for the interesting range
    for n in range(2, maxInteractions+2):
        pdfs[n] = getPdfOfSum(pdfs[n-1], pdfs[1], n) # add one interaction to previous pdf
        pdfs[n].Scale(1/pdfs[n].Integral()) # normalize
        pdfs[n].SetLineColor(n)
        if draw:
            pdfs[n].Draw("HISTSAME")
            ROOT.gPad.SetTitle("Track probabilities for different number of interactions;Track multiplicity;Probability")
    if wait:
        ROOT.gPad.WaitPrimitive()
    return pdfs

# thin wrapper for generating a function for the Poisson distribution describing the number of interactions for a given mu value
def getInteractionsPerBunchCrossingPdf(muNtrkPdfs, mu, muMax):#, wait = False, c = None):
    myPoisson = ROOT.TF1("myPois %f" % mu, "TMath::Poisson(x,[mu])", -0.5, muMax+0.5)
    myPoisson.SetParameter("mu", mu)
    myPoisson.SetNpx(muMax+1)
    return myPoisson

# we need to adopt a few limitations
measurementsMax = 200 # maximum number of tracks per BC we'll consider here
muMax = 80 # needs to be set in harmony with above, where the scaling depends on type of measurements used (and their typical number per interaction)

# now try it out - make a function
f = ROOT.TF1("myFunc1", "e**(-0.272581-1.933719*x)+e**(-2.107228-0.198227*x)", -0.5, 19.5)

c = ROOT.TCanvas("Lumi PDFs", "Lumi PDFs", 1200, 800)
c.Divide(2,2)
c.cd(1)

# get the per-interaction histogram (can also be skipped if already have a histo and don't need a function!)
c.Print("LumiPlots.pdf[")
h = getPdfHistogramFromFunction(f, True, c)
c.Print("LumiPlots.pdf")

c.cd(2)
pdfDict = makeNInteractionHistos(h, muMax, measurementsMax, True, c)
c.Print("LumiPlots.pdf")

# now loop over a couple of average mu values and construct the
nIntPerBcDict = {}
nTrkPerBcDict = {}
c.cd(3)
leg = ROOT.TLegend(0.72,0.63,0.88,0.85)
leg.SetBorderSize(0)
leg.SetFillStyle(0)

print("Will now generate distributions for three test mu values")
for i, mu in enumerate([24, 27, 30]):
    print("Average mu = %f" % mu)
    c.cd(3)
    nIntPerBcDict[mu] = getInteractionsPerBunchCrossingPdf(pdfDict, mu, muMax)
    drawOpt = ""
    if i > 0:
        drawOpt += "SAME"
    #nIntPerBcDict[mu].SetLineColor(i+1)
    ph = nIntPerBcDict[mu].GetHistogram()
    ph.SetLineColor(i+1)
    ph.Draw(drawOpt)
    leg.AddEntry(ph, "#mu = %d" % mu, "l")
    leg.Draw()
    if wait:
        ROOT.gPad.WaitPrimitive()

    # now weight together the distributions for fixed nInt to this average mu and draw the final plot, the expected nTrk per BC
    c.cd(4)
    nTrkPerBcDict[mu] = ROOT.TH1D("Tracks per bunch crossing for mu of %f" % mu, "Tracks per bunch crossing for <#mu>", measurementsMax+1, -0.5, measurementsMax+0.5)
    for nInt in range(0, ph.GetNbinsX()+1):
        nIntProb = ph.GetBinContent(nInt+1)
        #print("  Adding nTrk distribution for %d interactions with weight %f" % (nInt, nIntProb))
        nTrkPerBcDict[mu].Add(nIntProb*pdfDict[nInt])
    nTrkPerBcDict[mu].SetLineColor(i+1)
    nTrkPerBcDict[mu].Draw("HIST "+drawOpt)
    leg.Draw()
    c.Print("LumiPlots.pdf")

# define a function that can be used to fit to a measured distribution to determine mu
def tracksPerBC(x, params):
    x = x[0]
    mu = params[0]
    norm = params[1]
    #print("Will now evaluate tracksPerBC at x = %f with mu = %f" % (x, mu))
    # make a tracks-per-BC histogram for the given mu value
    # 1. get the Poisson distribution
    nIntPerBcFunction = getInteractionsPerBunchCrossingPdf(pdfDict, mu, muMax)
    #print(nIntPerBcFunction)
    nIntPerBcHisto = nIntPerBcFunction.GetHistogram()
    #print(nIntPerBcHisto)
    # 2. weight together the tracks-per-BC PDFs for fixed nInt
    nTrkPerBcHisto = ROOT.TH1D("TrackPDF", "TrackPDF", measurementsMax+1, -0.5, measurementsMax+0.5)
    for nInt in range(0, ph.GetNbinsX()+1):
        nIntProb = nIntPerBcHisto.GetBinContent(nInt+1)
        #print("%d: %f" % (nInt, nIntProb))
        nTrkPerBcHisto.Add(nIntProb*pdfDict[nInt])
    nTrkPerBcHisto.Scale(1/nTrkPerBcHisto.Integral())
    # now scale the PDF with the norm parameter
    nTrkPerBcHisto.Scale(norm)
    # return the value for the specified x (0 tracks is in bin 1)
    #print(nTrkPerBcHisto)
    return nTrkPerBcHisto.GetBinContent(int(x)+1)

# now define the TF1 object
tracksPerBCFunction = ROOT.TF1("TracksPerBC", tracksPerBC, -0.5, measurementsMax+0.5, 2) # last argument is the number of parameters of the function

# generate some data
mu = 12.4
nIntPerBc = getInteractionsPerBunchCrossingPdf(pdfDict, mu, muMax)
nIntPerBcHisto = nIntPerBc.GetHistogram()
nTracksPerBCHisto = ROOT.TH1D("Tracks per BC for mu of %f" % mu, "Tracks per bunch crossing for <#mu> = %f" % mu, measurementsMax+1, -0.5, measurementsMax+0.5)
for event in range(0, 10000):
    tracks = 0
    nInt = nIntPerBc.GetRandom()
    #print(nInt)
    for i in range(0, int(nInt)+1):
        tracks += pdfDict[1].GetRandom()
    nTracksPerBCHisto.Fill(tracks)
nTracksPerBCHisto.Draw()
print("Generated mock data, will now try to fit it")

# now we need to set good starting guesses for the fit parameters
# for the mean, we should start with the current way to determine the luminosity, i.e. just taking the sum of tracks and dividing by the mean number of tracks per interaction
guessedMu = nTracksPerBCHisto.GetMean()/pdfDict[1].GetMean()
guessedNorm = nTracksPerBCHisto.Integral()
tracksPerBCFunction.SetParameter(0, guessedMu)
tracksPerBCFunction.SetParameter(1, guessedNorm)
print("Will start with guesses: mu = %f, norm = %f" % (guessedMu, guessedNorm))
fr = nTracksPerBCHisto.Fit(tracksPerBCFunction, "S")
print("Chi2: %f, nDoF: %d --> Chi2/nDoF = %f" % (fr.Chi2(), fr.Ndf(), fr.Chi2()/fr.Ndf()))
c.Update()
fr.Print()

c.Print("LumiPlots.pdf")
c.Print("LumiPlots.pdf]") # close the file
