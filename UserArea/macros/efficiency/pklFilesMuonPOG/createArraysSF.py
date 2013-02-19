#!/usr/bin/env python

import pickle

# open the files
f_2012AB = open('MuonEfficiencies_Run_2012A_2012_B_53X.pkl', 'r')
f_2012C  = open('MuonEfficiencies_Run_2012C_53X.pkl', 'r')

# dictionary
dict_2012AB = pickle.load( f_2012AB)
dict_2012C  = pickle.load( f_2012C )


# defines the eta/pt ranges
etaRange = ['-2.1_-1.6',
            '-1.6_-1.2',
            '-1.2_-0.9',
            '-0.9_-0.6',
            '-0.6_-0.3',
            '-0.3_-0.2',
            '-0.2_0.2',
            '0.2_0.3',
            '0.3_0.6',
            '0.6_0.9',
            '0.9_1.2',                  
            '1.2_1.6',
            '1.6_2.1']

ptRange = ['25_30',
           '30_35',
           '35_40',
           '40_50',
           '50_60',
           '60_90',
           '90_140',
           '140_500']


print "const int etaBins = 13;"
print "const int ptBins = 8;"
print " "
print "double etaRange[etaBins+1] = {-2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2,"
print "                              0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1};"
print " "
print "double ptRange[ptBins+1] = {25, 30, 35, 40, 50, 60, 90, 140, 500};"
print " "
print " "

print "// 2012AB ID/Iso" 
print "double MuonIdTight_2012AB[etaBins] = {"
for eta in etaRange:
    print dict_2012AB['Tight']['etapt20-500'][eta]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double PFIsoIddB_eta09_2012AB[ptBins] = {"
for pt in ptRange:
    print dict_2012AB['combRelIsoPF04dBeta<012_Tight']['ptabseta<0.9'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double PFIsoIddB_eta09to12_2012AB[ptBins] = {"
for pt in ptRange:
    print dict_2012AB['combRelIsoPF04dBeta<012_Tight']['ptabseta0.9-1.2'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double PFIsoIddB_eta12to21_2012AB[ptBins] = {"
for pt in ptRange:
    print dict_2012AB['combRelIsoPF04dBeta<012_Tight']['ptabseta1.2-2.1'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "// Trigger efficiencies split by run periods"
print "double TrigId_2012A[etaBins] = {"
for eta in etaRange:
    print dict_2012AB['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500_2012A'][eta]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"

print "double TrigId_2012B[etaBins] = {"
for eta in etaRange:
    print dict_2012AB['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500_2012B'][eta]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"



print "// 2012C ID/Iso/Trigger" 
print "double MuonIdTight_2012C[etaBins] = {"
for eta in etaRange:
    print dict_2012C['Tight']['etapt20-500'][eta]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double PFIsoIddB_eta09_2012C[ptBins] = {"
for pt in ptRange:
    print dict_2012C['combRelIsoPF04dBeta<012_Tight']['ptabseta<0.9'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double PFIsoIddB_eta09to12_2012C[ptBins] = {"
for pt in ptRange:
    print dict_2012C['combRelIsoPF04dBeta<012_Tight']['ptabseta0.9-1.2'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double PFIsoIddB_eta12to21_2012C[ptBins] = {"
for pt in ptRange:
    print dict_2012C['combRelIsoPF04dBeta<012_Tight']['ptabseta1.2-2.1'][pt]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"


print "double TrigId_2012C[etaBins] = {"
for eta in etaRange:
    print dict_2012C['IsoMu24_eta2p1_TightIso']['eta_2p1pt25-500'][eta]['data/mc']['efficiency_ratio'], ", "
print "};\n\n"

