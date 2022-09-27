#!/usr/bin/env python3

from __future__ import print_function

import fastjet as fj
import fjcontrib
import fjext

import tqdm
import argparse
import os
import numpy as np
import array 
import math

from heppy.pythiautils import configuration as pyconf
import pythia8
import pythiafjext
import pythiaext

import ROOT

def calculate_distance(p0, p1):
    
	dphiabs = math.fabs(p0.phi() - p1.phi())
	dphi = dphiabs

	if dphiabs > math.pi:
		dphi = 2*math.pi - dphiabs


	dy = p0.rap() - p1.rap() # FIX ME: use p.eta() instead of p.rap()
	return math.sqrt(dy*dy + dphi*dphi)

class CorrelatorBuilder:
	def __init__(self, particle_list, scale):
		self.particle_list = particle_list
		self.mult = len(self.particle_list)
		self.pair_list = np.empty((0,self.mult),float)
		self.scale = scale

	def make_pairs(self):
		for i,part_i in enumerate(self.particle_list):
			inner_list = np.array([])
			for j,part_j in enumerate(self.particle_list[i:]):
				
				dist = calculate_distance(part_i,part_j)
				# print(' pairing particle i = ', i, ' and j = ', j+i, ' with distance ', dist)
				inner_list = np.append(inner_list, dist)

			inner_list = np.pad(inner_list, (i, 0), 'constant')
			inner_list.resize(1,self.mult)
			self.pair_list = np.append(self.pair_list,inner_list,axis=0)

			del inner_list

		# print(self.pair_list)

	def construct_EEC(self, hist):
		for ipart1 in range(self.mult):
			overlap = 0
			for ipart2 in range(ipart1,self.mult):
				dist12 = self.pair_list[ipart1][ipart2]
				# print(' EEC combining particle i =', ipart1, 'and j =', ipart2, 'with distance', dist12)
				
				# print(' E(i) =', self.particle_list[ipart1].E(), 'and E(j) =', self.particle_list[ipart2].E(), 'with jet pt =', self.scale)
				eec_weight = self.particle_list[ipart1].E()*self.particle_list[ipart2].E()/math.pow(self.scale,2) # NB: start with pt() instead of E()
				if ipart1 == ipart2:
					overlap = overlap + 1

				if overlap == 0:
				    eec_weight = 2*eec_weight
				if overlap > 0:
				    eec_weight = 1*eec_weight

				hist.Fill(dist12,eec_weight)

	def construct_E3C(self, hist):
		for ipart1 in range(self.mult):
			overlap = 0
			for ipart2 in range(ipart1,self.mult):
				dist12 = self.pair_list[ipart1][ipart2]
				if ipart1 == ipart2:
					overlap = overlap + 1
				for ipart3 in range(ipart2,self.mult):
					dist23 = self.pair_list[ipart2][ipart3]
					dist13 = self.pair_list[ipart1][ipart3]

					dist_list= [dist12, dist23, dist13]
					dist_list_sorted = sorted(dist_list)
					dist_max = dist_list_sorted[len(dist_list)-1]
					# print(' E3C combining particle', ipart1, ipart2, ipart3, 'with distance', dist12, dist23, dist13,'max',dist_max)
					
					if ipart2 == ipart3:
					    overlap = overlap + 1

					if ipart1 == ipart3:
					    overlap = overlap + 1

					e3c_weight = self.particle_list[ipart1].E()*self.particle_list[ipart2].E()*self.particle_list[ipart3].E()/math.pow(self.scale,3)

					if overlap == 0:
						e3c_weight = 3*2*e3c_weight
					if overlap == 1:
						e3c_weight = 3*e3c_weight
					if overlap > 1:
						e3c_weight = 1*e3c_weight

					hist.Fill(dist_max,e3c_weight)

	def construct_E4C(self, hist):
		for ipart1 in range(self.mult):
			overlap = 0
			for ipart2 in range(ipart1,self.mult):
				dist12 = self.pair_list[ipart1][ipart2]
				if ipart1 == ipart2:
					overlap = overlap + 1
				for ipart3 in range(ipart2,self.mult):
					dist23 = self.pair_list[ipart2][ipart3]
					dist13 = self.pair_list[ipart1][ipart3]
					if ipart2 == ipart3:
					    overlap = overlap + 1
					if ipart1 == ipart3:
					    overlap = overlap + 1

					for ipart4 in range(ipart3,self.mult):
						dist34 = self.pair_list[ipart3][ipart4]
						dist24 = self.pair_list[ipart2][ipart4]
						dist14 = self.pair_list[ipart1][ipart4]
						if ipart3 == ipart4:
							overlap = overlap + 1
						if ipart2 == ipart4:
							overlap = overlap + 1
						if ipart1 == ipart4:y
							overlap = overlap + 1

						dist_list= [dist12, dist23, dist13, dist34, dist24, dist14]
						dist_list_sorted = sorted(dist_list)
						dist_max = dist_list_sorted[len(dist_list)-1]
				
						e4c_weight = self.particle_list[ipart1].E()*self.particle_list[ipart2].E()*self.particle_list[ipart3].E()*self.particle_list[ipart4].E()/math.pow(self.scale,4)

						if overlap == 0:
							e4c_weight = 4*3*2*e4c_weight
						if overlap == 1:
							e4c_weight = 4*3*e4c_weight
						if overlap == 3:
							e4c_weight = 4*e4c_weight
						if overlap > 3:
							e4c_weight = 1*e4c_weight

						hist.Fill(dist_max,e4c_weight)

	def construct_E5C(self, hist):
		for ipart1 in range(self.mult):

			for ipart2 in range(ipart1+1,self.mult):
				dist12 = self.pair_list[ipart1][ipart2]

				for ipart3 in range(ipart2+1,self.mult):
					dist23 = self.pair_list[ipart2][ipart3]
					dist13 = self.pair_list[ipart1][ipart3]

					for ipart4 in range(ipart3+1,self.mult):
						dist34 = self.pair_list[ipart3][ipart4]
						dist24 = self.pair_list[ipart2][ipart4]
						dist14 = self.pair_list[ipart1][ipart4]

						for ipart5 in range(ipart4+1,self.mult):
							dist45 = self.pair_list[ipart4][ipart5]
							dist35 = self.pair_list[ipart3][ipart5]
							dist25 = self.pair_list[ipart2][ipart5]
							dist15 = self.pair_list[ipart1][ipart5]

							dist_list= [dist12, dist23, dist13, dist34, dist24, dist14, dist45, dist35, dist25, dist15]
							dist_list_sorted = sorted(dist_list)
							dist_max = dist_list_sorted[len(dist_list)-1]
					
							e5c_weight = self.particle_list[ipart1].E()*self.particle_list[ipart2].E()*self.particle_list[ipart3].E()*self.particle_list[ipart4].E()*self.particle_list[ipart5].E()/math.pow(self.scale,5)
							hist.Fill(dist_max,e5c_weight)


def linbins(xmin, xmax, nbins):
	lspace = np.linspace(xmin, xmax, nbins+1)
	arr = array.array('f', lspace)
	return arr

def logbins(xmin, xmax, nbins):
	lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
	arr = array.array('f', lspace)
	return arr

def main():
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly', prog=os.path.basename(__file__))
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--jet-R', help="jet radius", default=0.4, type=float)
	parser.add_argument('--min-jet-pt', help="minimum jet pT to accept", default=10.0, type=float)
	parser.add_argument('--max-jet-pt', help="maximum jet pT to accept", default=50.0, type=float)
	parser.add_argument('--max-jet-eta', help="maximum jet eta to accept", default=3.0, type=float)

	args = parser.parse_args()
	print('running with args:', args)

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_def = fj.JetDefinition(fj.antikt_algorithm, args.jet_R)
	jet_selector = fj.SelectorPtMin(args.min_jet_pt) & fj.SelectorAbsEtaMax(args.max_jet_eta)
	print(jet_def)

	all_jets = []

	mycfg = [] # no hardcoded settings
	pythia = pyconf.create_and_init_pythia_from_args(args, mycfg)

	link_fastjetvectors_to_pythia_particles = False
	user_index_start = 0

	if not pythia:
		print("[e] pythia initialization failed.")
		return
	if args.nev < 2:
		args.nev = 2
        
    # Create a histogram with ROOT
	lbins1 = linbins(0., 800, 800)
	hJetPt = ROOT.TH1D("hJetPt", "hJetPt", 800, lbins1)
	lbins2 = logbins(1E-4,1,50)
	hJetEEC = ROOT.TH1D("hJetEEC", "hJetEEC", 50, lbins2)
	hJetE3C = ROOT.TH1D("hJetE3C", "hJetE3C", 50, lbins2)
	hJetE4C = ROOT.TH1D("hJetE4C", "hJetE4C", 50, lbins2)
	# hJetE5C = ROOT.TH1D("hJetE5C", "hJetE5C", 50, lbins2)
	
	# event loop
	for i in tqdm.tqdm(range(args.nev)):
		if not pythia.next():
			continue

		# transform pythia particles into fastjet::PseudoJet(s) - use only final state (after hadronization in this case)
		parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kFinal], user_index_start, link_fastjetvectors_to_pythia_particles)

		# parts = pythiafjext.vectorize_select(pythia, [pythiafjext.kCharged], user_index_start, link_fastjetvectors_to_pythia_particles)
		parts_selected = []
		for part in parts:
			if part.pt() > 0 and part.eta() < 2:
				parts_selected.append(part)
		jets = fj.sorted_by_pt(jet_selector(jet_def(parts_selected)))

		# loop over jets
		for j in jets:
			if j.perp() < 500:
				continue
			if j.perp() > 550:
				continue
			if math.fabs(j.eta()) > 1.9:
				continue
			print('jet [px, py, pz, E]:', j)
			print('    [pt, eta, phi] :', [j.perp(), j.eta(), j.phi()])
			hJetPt.Fill(j.perp())

			constituents = fj.sorted_by_pt(j.constituents())

			part_pt_thrd = 1.0
			c_select = []
			for c in constituents:
				if c.pt() > part_pt_thrd:
					ip = c.user_index() # find the index in the pythia particle list
					if pythia.event[ip].charge()!=0:
						# print('part index',ip,'px',c.px(),'vs pythia px',pythia.event[ip].px(),'charge',pythia.event[ip].charge())
						c_select.append(c)

					c_select.append(c)
					
					# print('part pt',c.pt())

			new_corr = CorrelatorBuilder(c_select,j.perp())
			new_corr.make_pairs()
			new_corr.construct_EEC(hJetEEC)
			new_corr.construct_E3C(hJetE3C)
			new_corr.construct_E4C(hJetE4C)
			# new_corr.construct_E5C(hJetE5C)

		print('end of the event processing...')

	pythia.stat()
    
	fout = ROOT.TFile("test.root", "recreate")
	fout.cd()
	hJetPt.Write()
	hJetEEC.Write()
	hJetE3C.Write()
	hJetE4C.Write()
	# hJetE5C.Write()
	fout.Close()

if __name__ == '__main__':
	main()
