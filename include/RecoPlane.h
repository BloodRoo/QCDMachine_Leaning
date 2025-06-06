#include "JetBranch.h"
#include "Matching.cc"
#include "ParticleInfo.h"
#include "Process_bar.h"
#include "TLorentzVector.h"
#include "TreeEvents.h"
#include "common_tool.h"
#include "fastjet/contrib/IFNPlugin.hh"
#include "fastjet/contrib/GHSAlgo.hh"
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <algorithm>
#include <cassert>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <variant>
#include <vector>
#include "JetObservable.h"
#include "SpinObservable.h"
using namespace fastjet;
using namespace std;
using namespace fastjet::contrib;

class RecoPlane
{
public:
  static JetBranch::twoplanes
  JetConstituents(vector<ParticleInfo> constituents)
  {
    double z1cut = 0;
    double z2cut = 0;
    double kt1cut = 1;
    double kt2cut = 0;
    double issecondsoft = true;
    vector<JetBranch::twoplanes> twoplanes;
    std::vector<PseudoJet> particles;
    std::vector<ParticleInfo> particlesinfo;
    for (size_t i = 0; i < constituents.size(); ++i)
    {
      TLorentzVector p;
      p.SetPtEtaPhiE(constituents.at(i).pt, constituents.at(i).eta,
                     constituents.at(i).phi, constituents.at(i).e);
      PseudoJet particle = PseudoJet(p.Px(), p.Py(), p.Pz(), p.Energy());
      int pdgid = constituents.at(i).pdgid;
      ParticleInfo particleInfo(pdgid, constituents.at(i).chargeInt);
      particlesinfo.push_back(particleInfo);
      particle.set_user_index(i);
      particle.set_user_info(new FlavHistory(pdgid));
      particles.push_back(particle);
    }
    JetDefinition jet_def_CA(cambridge_algorithm,
                             JetDefinition::max_allowable_R);
    fastjet::contrib::FlavRecombiner flav_recombiner;
    jet_def_CA.set_recombiner(&flav_recombiner);
    double alpha = 1.0;
    double omega = 3.0 - alpha;
    fastjet::contrib::FlavRecombiner::FlavSummation flav_summation =
        fastjet::contrib::FlavRecombiner::net;
    auto ifn_plugin = new IFNPlugin(jet_def_CA, alpha, omega, flav_summation);
    JetDefinition IFN_jet_def(ifn_plugin);
    IFN_jet_def.delete_plugin_when_unused();
    ClusterSequence cs_IFN(particles, IFN_jet_def);
    vector<PseudoJet> jets_IFN = sorted_by_pt(cs_IFN.inclusive_jets());
    PseudoJet j0, j1, j2, j3, j4;
    j0 = jets_IFN[0];

    double Q = j0.e();
    auto twoplane = JetBranch::findPrimaryAndSecondaryJets(
        j0, z1cut, z2cut, kt1cut * 2 * Q, kt2cut * 2 * Q, particlesinfo,
        issecondsoft);

    return twoplane;
  }
  static JetBranch::threeplanes
  JetConstituents_threeplanes(vector<ParticleInfo> constituents)
  {
    double z1cut = 0;
    double z2cut = 0;
    double kt1cut = 1;
    double kt2cut = 0;
    double issecondsoft = true;
    vector<JetBranch::twoplanes> twoplanes;
    std::vector<PseudoJet> particles;
    std::vector<ParticleInfo> particlesinfo;
    int bnum = 0;
    for (size_t i = 0; i < constituents.size(); ++i)
    {
      TLorentzVector p;
      p.SetPtEtaPhiE(constituents.at(i).pt, constituents.at(i).eta,
                     constituents.at(i).phi, constituents.at(i).e);
      PseudoJet particle = PseudoJet(p.Px(), p.Py(), p.Pz(), p.Energy());
      int pdgid = constituents.at(i).pdgid;
      ParticleInfo particleInfo(pdgid, constituents.at(i).chargeInt);
      particlesinfo.push_back(particleInfo);
      particle.set_user_index(i);
      particle.set_user_info(new FlavHistory(pdgid));
      particles.push_back(particle);
      if ((abs(pdgid / 100 % 10) == 5 ||
           abs(pdgid / 1000 % 10) == 5))
        bnum++;
    }
    JetDefinition jet_def_CA(cambridge_algorithm,
                             JetDefinition::max_allowable_R);
    fastjet::contrib::FlavRecombiner flav_recombiner;
    jet_def_CA.set_recombiner(&flav_recombiner);

    double alpha = 1.0;
    double omega = 3.0 - alpha;
    fastjet::contrib::FlavRecombiner::FlavSummation flav_summation =
        fastjet::contrib::FlavRecombiner::net;
    auto ifn_plugin = new IFNPlugin(jet_def_CA, alpha, omega, flav_summation);
    JetDefinition IFN_jet_def(ifn_plugin);
    IFN_jet_def.delete_plugin_when_unused();
    ClusterSequence cs_IFN(particles, IFN_jet_def);
    vector<PseudoJet> jets_IFN = sorted_by_pt(cs_IFN.inclusive_jets());
    PseudoJet j0, j1, j2, j3, j4;
    j0 = jets_IFN[0];
    double Q = j0.e();
    auto threeplanes = JetBranch::findPrimarySecondaryAndThirdaryJets(
        j0, z1cut, z2cut, kt1cut * 2 * Q, kt2cut * 2 * Q, particlesinfo,
        issecondsoft);
    return threeplanes;
  }

  static fastjet::PseudoJet CAJet(vector<ParticleInfo> constituents)
  {
    std::vector<PseudoJet> particles;
    std::vector<ParticleInfo> particlesinfo;
    int bnum = 0;
    for (size_t i = 0; i < constituents.size(); ++i)
    {
      TLorentzVector p;
      p.SetPtEtaPhiE(constituents.at(i).pt, constituents.at(i).eta,
                     constituents.at(i).phi, constituents.at(i).e);
      PseudoJet particle = PseudoJet(p.Px(), p.Py(), p.Pz(), p.Energy());
      int pdgid = constituents.at(i).pdgid;
      ParticleInfo particleInfo(pdgid, constituents.at(i).chargeInt);
      particlesinfo.push_back(particleInfo);
      particle.set_user_index(i);
      particle.set_user_info(new FlavHistory(pdgid));
      particles.push_back(particle);
      if ((abs(pdgid / 100 % 10) == 5 ||
           abs(pdgid / 1000 % 10) == 5))
        bnum++;
    }
    JetDefinition jet_def_CA(cambridge_algorithm,
                             JetDefinition::max_allowable_R);
    fastjet::contrib::FlavRecombiner flav_recombiner;
    jet_def_CA.set_recombiner(&flav_recombiner);

    double alpha = 1.0;
    double omega = 3.0 - alpha;
    fastjet::contrib::FlavRecombiner::FlavSummation flav_summation =
        fastjet::contrib::FlavRecombiner::net;
    auto ifn_plugin = new IFNPlugin(jet_def_CA, alpha, omega, flav_summation);
    JetDefinition IFN_jet_def(ifn_plugin);
    IFN_jet_def.delete_plugin_when_unused();
    ClusterSequence cs_IFN(particles, IFN_jet_def);
    vector<PseudoJet> jets_IFN = sorted_by_pt(cs_IFN.inclusive_jets());
    PseudoJet j0, j1, j2, j3, j4;
    j0 = jets_IFN[0];

    return j0;
  }
  // static JetBranch::BhadronPlanes
  // JetConstituents_bhadronplanes(vector<ParticleInfo> constituents)
  // {
  //   double z1cut = 0;
  //   double z2cut = 0;
  //   double kt1cut = 1;
  //   double kt2cut = 0;
  //   double issecondsoft = true;
  //   vector<JetBranch::twoplanes> twoplanes;
  //   std::vector<PseudoJet> particles;
  //   std::vector<ParticleInfo> particlesinfo;
  //   for (size_t i = 0; i < constituents.size(); ++i)
  //   {
  //     TLorentzVector p;
  //     p.SetPtEtaPhiE(constituents.at(i).pt, constituents.at(i).eta,
  //                    constituents.at(i).phi, constituents.at(i).e);
  //     PseudoJet particle = PseudoJet(p.Px(), p.Py(), p.Pz(), p.Energy());
  //     int pdgid = constituents.at(i).pdgid;
  //     ParticleInfo particleInfo(pdgid, constituents.at(i).chargeInt);
  //     particlesinfo.push_back(particleInfo);
  //     particle.set_user_index(i);
  //     particle.set_user_info(new FlavHistory(pdgid));
  //     particles.push_back(particle);
  //   }
  //   JetDefinition jet_def_CA(cambridge_algorithm,
  //                            JetDefinition::max_allowable_R);
  //   fastjet::contrib::FlavRecombiner flav_recombiner;
  //   jet_def_CA.set_recombiner(&flav_recombiner);

  //   double alpha = 1.0;
  //   double omega = 3.0 - alpha;
  //   fastjet::contrib::FlavRecombiner::FlavSummation flav_summation =
  //       fastjet::contrib::FlavRecombiner::net;
  //   auto ifn_plugin = new IFNPlugin(jet_def_CA, alpha, omega, flav_summation);
  //   JetDefinition IFN_jet_def(ifn_plugin);
  //   IFN_jet_def.delete_plugin_when_unused();
  //   ClusterSequence cs_IFN(particles, IFN_jet_def);
  //   vector<PseudoJet> jets_IFN = sorted_by_pt(cs_IFN.inclusive_jets());
  //   PseudoJet j0, j1, j2, j3, j4;
  //   j0 = jets_IFN[0];
  //   double Q = j0.e();
  //   auto threeplanes = JetBranch::findBhadronPlanes(j0, particlesinfo);

  //   return threeplanes;
  // }
  static JetBranch::threeplanes
  JetConstituents_threeplanesGHS(vector<ParticleInfo> constituents)
  {
    double z1cut = 0;
    double z2cut = 0;
    double kt1cut = 1;
    double kt2cut = 0;
    double issecondsoft = true;
    vector<JetBranch::twoplanes> twoplanes;
    std::vector<PseudoJet> particles;
    std::vector<ParticleInfo> particlesinfo;
    for (size_t i = 0; i < constituents.size(); ++i)
    {
      TLorentzVector p;
      p.SetPtEtaPhiE(constituents.at(i).pt, constituents.at(i).eta,
                     constituents.at(i).phi, constituents.at(i).e);
      PseudoJet particle = PseudoJet(p.Px(), p.Py(), p.Pz(), p.Energy());
      int pdgid = constituents.at(i).pdgid;
      ParticleInfo particleInfo(pdgid, constituents.at(i).chargeInt);
      particlesinfo.push_back(particleInfo);
      particle.set_user_index(i);
      particle.set_user_info(new FlavHistory(pdgid));
      particles.push_back(particle);
    }
    JetDefinition jet_def_CA(cambridge_algorithm,
                             JetDefinition::max_allowable_R);
    fastjet::contrib::FlavRecombiner flav_recombiner;
    jet_def_CA.set_recombiner(&flav_recombiner);

    double GHS_alpha = 1.0; // < flav-kt distance parameter alpha
    double GHS_omega = 2.0;
    ClusterSequence cs_CA(particles, jet_def_CA);
    vector<PseudoJet> jets_CA = sorted_by_pt(cs_CA.inclusive_jets());
    vector<PseudoJet> GHS_jets = run_GHS(jets_CA, 0, GHS_alpha, GHS_omega, flav_recombiner);

    PseudoJet j0, j1, j2, j3, j4;
    j0 = GHS_jets[0];

    double Q = j0.e();
    auto threeplanes = JetBranch::findPrimarySecondaryAndThirdaryJets(
        j0, z1cut, z2cut, kt1cut * 2 * Q, kt2cut * 2 * Q, particlesinfo,
        issecondsoft);

    return threeplanes;
  }
  static void SavePlanes(JetBranch::twoplanes twoplanes, TreeEvents &treeEvents,
                         int jetindex, string prefix = "", string suffix = "")
  {
    treeEvents.push_back(prefix + "pt1" + suffix, twoplanes.first.harder.pt());
    treeEvents.push_back(prefix + "eta1" + suffix,
                         twoplanes.first.harder.eta());
    treeEvents.push_back(prefix + "phi1" + suffix,
                         twoplanes.first.harder.phi());
    treeEvents.push_back(prefix + "e1" + suffix, twoplanes.first.harder.e());
    treeEvents.push_back(prefix + "ntracks1" + suffix,
                         twoplanes.first.harder_ntracks);
    treeEvents.push_back(prefix + "nparticles1" + suffix,
                         twoplanes.first.harder_nparticles);

    treeEvents.push_back(prefix + "pt2" + suffix, twoplanes.first.softer.pt());
    treeEvents.push_back(prefix + "eta2" + suffix,
                         twoplanes.first.softer.eta());
    treeEvents.push_back(prefix + "phi2" + suffix,
                         twoplanes.first.softer.phi());
    treeEvents.push_back(prefix + "e2" + suffix, twoplanes.first.softer.e());
    treeEvents.push_back(prefix + "ntracks2" + suffix,
                         twoplanes.first.softer_ntracks);
    treeEvents.push_back(prefix + "nparticles2" + suffix,
                         twoplanes.first.softer_nparticles);

    treeEvents.push_back(prefix + "pt3" + suffix, twoplanes.second.harder.pt());
    treeEvents.push_back(prefix + "eta3" + suffix,
                         twoplanes.second.harder.eta());
    treeEvents.push_back(prefix + "phi3" + suffix,
                         twoplanes.second.harder.phi());
    treeEvents.push_back(prefix + "e3" + suffix, twoplanes.second.harder.e());
    treeEvents.push_back(prefix + "ntracks3" + suffix,
                         twoplanes.second.harder_ntracks);
    treeEvents.push_back(prefix + "nparticles3" + suffix,
                         twoplanes.second.harder_nparticles);

    treeEvents.push_back(prefix + "pt4" + suffix, twoplanes.second.softer.pt());
    treeEvents.push_back(prefix + "eta4" + suffix,
                         twoplanes.second.softer.eta());
    treeEvents.push_back(prefix + "phi4" + suffix,
                         twoplanes.second.softer.phi());
    treeEvents.push_back(prefix + "e4" + suffix, twoplanes.second.softer.e());
    treeEvents.push_back(prefix + "ntracks4" + suffix,
                         twoplanes.second.softer_ntracks);
    treeEvents.push_back(prefix + "nparticles4" + suffix,
                         twoplanes.second.softer_nparticles);

    treeEvents.push_back(prefix + "z1" + suffix, twoplanes.first.z);
    treeEvents.push_back(prefix + "kt1" + suffix, twoplanes.first.kt);
    treeEvents.push_back(prefix + "theta1" + suffix, twoplanes.first.theta);
    treeEvents.push_back(prefix + "deltaR1" + suffix, twoplanes.first.deltaR);
    treeEvents.push_back(
        prefix + "eec1" + suffix,
        (twoplanes.first.harder.e() * twoplanes.first.softer.e()) * 1.0 /
            twoplanes.Q / twoplanes.Q);
    treeEvents.push_back(prefix + "z2" + suffix, twoplanes.second.z);
    treeEvents.push_back(prefix + "kt2" + suffix, twoplanes.second.kt);
    treeEvents.push_back(prefix + "theta2" + suffix, twoplanes.second.theta);
    treeEvents.push_back(prefix + "deltaR2" + suffix, twoplanes.second.deltaR);
    treeEvents.push_back(
        prefix + "eec2" + suffix,
        (twoplanes.second.harder.e() * twoplanes.second.softer.e()) * 1.0 /
            twoplanes.Q / twoplanes.Q);
    treeEvents.push_back(prefix + "n1x" + suffix, twoplanes.first.n.x());
    treeEvents.push_back(prefix + "n1y" + suffix, twoplanes.first.n.y());
    treeEvents.push_back(prefix + "n1z" + suffix, twoplanes.first.n.z());
    treeEvents.push_back(prefix + "n2x" + suffix, twoplanes.second.n.x());
    treeEvents.push_back(prefix + "n2y" + suffix, twoplanes.second.n.y());
    treeEvents.push_back(prefix + "n2z" + suffix, twoplanes.second.n.z());
    treeEvents.push_back(prefix + "isqq" + suffix,
                         static_cast<int>(twoplanes.second.isqq));
    treeEvents.push_back(prefix + "isgg" + suffix,
                         static_cast<int>(twoplanes.second.isgg));
    treeEvents.push_back(prefix + "phi" + suffix, twoplanes.deltaPhi);

    treeEvents.push_back(prefix + "jpt" + suffix, twoplanes.first.initJet.Pt());
    treeEvents.push_back(prefix + "jeta" + suffix,
                         twoplanes.first.initJet.Eta());
    treeEvents.push_back(prefix + "jphi" + suffix,
                         twoplanes.first.initJet.Phi());
    treeEvents.push_back(prefix + "je" + suffix, twoplanes.first.initJet.E());

    for (int ii = 0; ii < twoplanes.first.harder_constituents_info.size();
         ii++)
    {
      auto particle = twoplanes.first.harder_constituents_info.at(ii);
      treeEvents.push_back(prefix + "particle1_pt" + suffix, particle.pt);
      treeEvents.push_back(prefix + "particle1_eta" + suffix, particle.eta);
      treeEvents.push_back(prefix + "particle1_phi" + suffix, particle.phi);
      treeEvents.push_back(prefix + "particle1_e" + suffix, particle.e);
      treeEvents.push_back(prefix + "particle1_charge" + suffix,
                           particle.charge);
      treeEvents.push_back(prefix + "particle1_pid" + suffix, particle.pdgid);
      treeEvents.push_back(prefix + "particle1_jetid" + suffix, jetindex);
    }
    for (int ii = 0; ii < twoplanes.first.softer_constituents_info.size();
         ii++)
    {
      auto particle = twoplanes.first.softer_constituents_info.at(ii);
      treeEvents.push_back(prefix + "particle2_pt" + suffix, particle.pt);
      treeEvents.push_back(prefix + "particle2_eta" + suffix, particle.eta);
      treeEvents.push_back(prefix + "particle2_phi" + suffix, particle.phi);
      treeEvents.push_back(prefix + "particle2_e" + suffix, particle.e);
      treeEvents.push_back(prefix + "particle2_charge" + suffix,
                           particle.charge);
      treeEvents.push_back(prefix + "particle2_pid" + suffix, particle.pdgid);
      treeEvents.push_back(prefix + "particle2_jetid" + suffix, jetindex);
    }

    for (int ii = 0; ii < twoplanes.second.harder_constituents_info.size();
         ii++)
    {
      auto particle = twoplanes.second.harder_constituents_info.at(ii);
      treeEvents.push_back(prefix + "particle3_pt" + suffix, particle.pt);
      treeEvents.push_back(prefix + "particle3_eta" + suffix, particle.eta);
      treeEvents.push_back(prefix + "particle3_phi" + suffix, particle.phi);
      treeEvents.push_back(prefix + "particle3_e" + suffix, particle.e);
      treeEvents.push_back(prefix + "particle3_charge" + suffix,
                           particle.charge);
      treeEvents.push_back(prefix + "particle3_pid" + suffix, particle.pdgid);
      treeEvents.push_back(prefix + "particle3_jetid" + suffix, jetindex);
    }
    for (int ii = 0; ii < twoplanes.second.softer_constituents_info.size();
         ii++)
    {
      auto particle = twoplanes.second.softer_constituents_info.at(ii);
      treeEvents.push_back(prefix + "particle4_pt" + suffix, particle.pt);
      treeEvents.push_back(prefix + "particle4_eta" + suffix, particle.eta);
      treeEvents.push_back(prefix + "particle4_phi" + suffix, particle.phi);
      treeEvents.push_back(prefix + "particle4_e" + suffix, particle.e);
      treeEvents.push_back(prefix + "particle4_charge" + suffix,
                           particle.charge);
      treeEvents.push_back(prefix + "particle4_pid" + suffix, particle.pdgid);
      treeEvents.push_back(prefix + "particle4_jetid" + suffix, jetindex);
    }
  }
  static void SaveLoopPlanes(JetBranch::threeplanes threeplanes, TreeEvents &treeEvents,
                             int jetindex, string prefix = "", string suffix = "",
                             bool SaveParticles = true, std::string SampleType = "PrivateMC")
  {
    treeEvents.push_back(prefix + "init_loop_pdgid3" + suffix,
                         threeplanes.second.harder_init_pdgid);
    treeEvents.push_back(prefix + "init_loop_descri3" + suffix,
                         threeplanes.second.harder_init_descri);
    treeEvents.push_back(prefix + "final_loop_pdgid3" + suffix,
                         threeplanes.second.harder_final_pdgid);
    treeEvents.push_back(prefix + "final_loop_descri3" + suffix,
                         threeplanes.second.harder_final_descri);

    treeEvents.push_back(prefix + "init_loop_pdgid4" + suffix,
                         threeplanes.second.softer_init_pdgid);
    treeEvents.push_back(prefix + "init_loop_descri4" + suffix,
                         threeplanes.second.softer_init_descri);
    treeEvents.push_back(prefix + "final_loop_pdgid4" + suffix,
                         threeplanes.second.softer_final_pdgid);
    treeEvents.push_back(prefix + "final_loop_descri4" + suffix,
                         threeplanes.second.softer_final_descri);

    treeEvents.push_back(prefix + "typeloop2" + suffix, JetObservable::determineType(threeplanes.second.softer_flav, threeplanes.second.harder_flav));
  }
  static void SavePlanes(JetBranch::threeplanes threeplanes, TreeEvents &treeEvents,
                         int jetindex, string prefix = "", string suffix = "",
                         bool SaveParticles = true, std::string SampleType = "PrivateMC")
  {
    JetObservable JetObservable1(threeplanes.first.harder_constituents_info);
    JetObservable JetObservable2(threeplanes.first.softer_constituents_info);
    JetObservable JetObservable3(threeplanes.second.harder_constituents_info);
    JetObservable JetObservable4(threeplanes.second.softer_constituents_info);
    JetObservable JetObservable5(threeplanes.third.harder_constituents_info);
    JetObservable JetObservable6(threeplanes.third.softer_constituents_info);
    auto width1 = JetObservable1.jetwidth();
    auto width2 = JetObservable2.jetwidth();
    auto width3 = JetObservable3.jetwidth();
    auto width4 = JetObservable4.jetwidth();
    auto width5 = JetObservable5.jetwidth();
    auto width6 = JetObservable6.jetwidth();

    JetObservable JetObservable31(threeplanes.second.harder_harder_constituents_info);
    JetObservable JetObservable32(threeplanes.second.harder_softer_constituents_info);
    JetObservable JetObservable41(threeplanes.second.softer_harder_constituents_info);
    JetObservable JetObservable42(threeplanes.second.softer_softer_constituents_info);
    JetObservable JetObservable51(threeplanes.third.harder_harder_constituents_info);
    JetObservable JetObservable52(threeplanes.third.harder_softer_constituents_info);
    JetObservable JetObservable61(threeplanes.third.softer_harder_constituents_info);
    JetObservable JetObservable62(threeplanes.third.softer_softer_constituents_info);
    auto width31 = JetObservable31.jetwidth();
    auto width32 = JetObservable32.jetwidth();
    auto width41 = JetObservable41.jetwidth();
    auto width42 = JetObservable42.jetwidth();
    auto width51 = JetObservable51.jetwidth();
    auto width52 = JetObservable52.jetwidth();
    auto width61 = JetObservable61.jetwidth();
    auto width62 = JetObservable62.jetwidth();

    SpinObservable spinobservable2(threeplanes.first.harder_constituents_info,
                                   threeplanes.first.softer_constituents_info,
                                   threeplanes.second.harder_constituents_info,
                                   threeplanes.second.softer_constituents_info);
    auto planetheta2 = spinobservable2.GetPlaneTheta();
    SpinObservable spinobservable3(threeplanes.first.softer_constituents_info,
                                   threeplanes.first.harder_constituents_info,
                                   threeplanes.third.harder_constituents_info,
                                   threeplanes.third.softer_constituents_info);
    auto planetheta3 = spinobservable3.GetPlaneTheta();

    treeEvents.push_back(prefix + "primaryindex" + suffix, threeplanes.firstindex);

    treeEvents.push_back(prefix + "pt1" + suffix, threeplanes.first.harder.pt());
    treeEvents.push_back(prefix + "eta1" + suffix,
                         threeplanes.first.harder.eta());
    treeEvents.push_back(prefix + "phi1" + suffix,
                         threeplanes.first.harder.phi());
    treeEvents.push_back(prefix + "e1" + suffix, threeplanes.first.harder.e());
    treeEvents.push_back(prefix + "ntracks1" + suffix,
                         threeplanes.first.harder_ntracks);
    treeEvents.push_back(prefix + "nparticles1" + suffix,
                         threeplanes.first.harder_nparticles);
    treeEvents.push_back(prefix + "flavour1" + suffix,
                         threeplanes.first.harder_flav);
    treeEvents.push_back(prefix + "pTD1" + suffix, JetObservable1.pTD());
    treeEvents.push_back(prefix + "sigma1" + suffix, width1.sigma);
    treeEvents.push_back(prefix + "sigma11" + suffix, width1.sigma1);
    treeEvents.push_back(prefix + "sigma21" + suffix, width1.sigma2);

    treeEvents.push_back(prefix + "pt2" + suffix, threeplanes.first.softer.pt());
    treeEvents.push_back(prefix + "eta2" + suffix,
                         threeplanes.first.softer.eta());
    treeEvents.push_back(prefix + "phi2" + suffix,
                         threeplanes.first.softer.phi());
    treeEvents.push_back(prefix + "e2" + suffix, threeplanes.first.softer.e());
    treeEvents.push_back(prefix + "ntracks2" + suffix,
                         threeplanes.first.softer_ntracks);
    treeEvents.push_back(prefix + "nparticles2" + suffix,
                         threeplanes.first.softer_nparticles);
    treeEvents.push_back(prefix + "flavour2" + suffix,
                         threeplanes.first.softer_flav);
    treeEvents.push_back(prefix + "pTD2" + suffix, JetObservable2.pTD());
    treeEvents.push_back(prefix + "sigma2" + suffix, width2.sigma);
    treeEvents.push_back(prefix + "sigma12" + suffix, width2.sigma1);
    treeEvents.push_back(prefix + "sigma22" + suffix, width2.sigma2);

    treeEvents.push_back(prefix + "pt3" + suffix, threeplanes.second.harder.pt());
    treeEvents.push_back(prefix + "eta3" + suffix,
                         threeplanes.second.harder.eta());
    treeEvents.push_back(prefix + "phi3" + suffix,
                         threeplanes.second.harder.phi());
    treeEvents.push_back(prefix + "e3" + suffix, threeplanes.second.harder.e());
    treeEvents.push_back(prefix + "ntracks3" + suffix,
                         threeplanes.second.harder_ntracks);
    treeEvents.push_back(prefix + "nparticles3" + suffix,
                         threeplanes.second.harder_nparticles);
    treeEvents.push_back(prefix + "flavour3" + suffix,
                         threeplanes.second.harder_flav);
    treeEvents.push_back(prefix + "pTD3" + suffix, JetObservable3.pTD());
    treeEvents.push_back(prefix + "sigma3" + suffix, width3.sigma);
    treeEvents.push_back(prefix + "sigma13" + suffix, width3.sigma1);
    treeEvents.push_back(prefix + "sigma23" + suffix, width3.sigma2);

    treeEvents.push_back(prefix + "pt4" + suffix, threeplanes.second.softer.pt());
    treeEvents.push_back(prefix + "eta4" + suffix,
                         threeplanes.second.softer.eta());
    treeEvents.push_back(prefix + "phi4" + suffix,
                         threeplanes.second.softer.phi());
    treeEvents.push_back(prefix + "e4" + suffix, threeplanes.second.softer.e());
    treeEvents.push_back(prefix + "ntracks4" + suffix,
                         threeplanes.second.softer_ntracks);
    treeEvents.push_back(prefix + "nparticles4" + suffix,
                         threeplanes.second.softer_nparticles);
    treeEvents.push_back(prefix + "flavour4" + suffix,
                         threeplanes.second.softer_flav);
    treeEvents.push_back(prefix + "pTD4" + suffix, JetObservable4.pTD());
    treeEvents.push_back(prefix + "sigma4" + suffix, width4.sigma);
    treeEvents.push_back(prefix + "sigma14" + suffix, width4.sigma1);
    treeEvents.push_back(prefix + "sigma24" + suffix, width4.sigma2);

    treeEvents.push_back(prefix + "pt5" + suffix, threeplanes.third.harder.pt());
    treeEvents.push_back(prefix + "eta5" + suffix,
                         threeplanes.third.harder.eta());
    treeEvents.push_back(prefix + "phi5" + suffix,
                         threeplanes.third.harder.phi());
    treeEvents.push_back(prefix + "e5" + suffix, threeplanes.third.harder.e());
    treeEvents.push_back(prefix + "ntracks5" + suffix,
                         threeplanes.third.harder_ntracks);
    treeEvents.push_back(prefix + "nparticles5" + suffix,
                         threeplanes.third.harder_nparticles);
    treeEvents.push_back(prefix + "flavour5" + suffix,
                         threeplanes.third.harder_flav);
    treeEvents.push_back(prefix + "pTD5" + suffix, JetObservable5.pTD());
    treeEvents.push_back(prefix + "sigma5" + suffix, width5.sigma);
    treeEvents.push_back(prefix + "sigma15" + suffix, width5.sigma1);
    treeEvents.push_back(prefix + "sigma25" + suffix, width5.sigma2);

    treeEvents.push_back(prefix + "pt6" + suffix, threeplanes.third.softer.pt());
    treeEvents.push_back(prefix + "eta6" + suffix,
                         threeplanes.third.softer.eta());
    treeEvents.push_back(prefix + "phi6" + suffix,
                         threeplanes.third.softer.phi());
    treeEvents.push_back(prefix + "e6" + suffix, threeplanes.third.softer.e());
    treeEvents.push_back(prefix + "ntracks6" + suffix,
                         threeplanes.third.softer_ntracks);
    treeEvents.push_back(prefix + "nparticles6" + suffix,
                         threeplanes.third.softer_nparticles);
    treeEvents.push_back(prefix + "flavour6" + suffix,
                         threeplanes.third.softer_flav);
    treeEvents.push_back(prefix + "pTD6" + suffix, JetObservable6.pTD());
    treeEvents.push_back(prefix + "sigma6" + suffix, width6.sigma);
    treeEvents.push_back(prefix + "sigma16" + suffix, width6.sigma1);
    treeEvents.push_back(prefix + "sigma26" + suffix, width6.sigma2);

    // treeEvents.push_back(prefix + "pt31" + suffix, threeplanes.second.harder_harder.pt());
    // treeEvents.push_back(prefix + "eta31" + suffix,
    //                      threeplanes.second.harder_harder.eta());
    // treeEvents.push_back(prefix + "phi31" + suffix,
    //                      threeplanes.second.harder_harder.phi());
    // treeEvents.push_back(prefix + "e31" + suffix, threeplanes.second.harder_harder.e());
    // treeEvents.push_back(prefix + "ntracks31" + suffix,
    //                      threeplanes.second.harder_harder_ntracks);
    // treeEvents.push_back(prefix + "nparticles31" + suffix,
    //                      threeplanes.second.harder_harder_nparticles);
    // treeEvents.push_back(prefix + "flavour31" + suffix,
    //                      JetBranch::GetIFNFlavour(threeplanes.second.harder_harder));
    // treeEvents.push_back(prefix + "pTD31" + suffix, JetObservable31.pTD());
    // treeEvents.push_back(prefix + "sigma31" + suffix, width31.sigma);
    // treeEvents.push_back(prefix + "sigma131" + suffix, width31.sigma1);
    // treeEvents.push_back(prefix + "sigma231" + suffix, width31.sigma2);

    // treeEvents.push_back(prefix + "pt32" + suffix, threeplanes.second.harder_softer.pt());
    // treeEvents.push_back(prefix + "eta32" + suffix,
    //                      threeplanes.second.harder_softer.eta());
    // treeEvents.push_back(prefix + "phi32" + suffix,
    //                      threeplanes.second.harder_softer.phi());
    // treeEvents.push_back(prefix + "e32" + suffix, threeplanes.second.harder_softer.e());
    // treeEvents.push_back(prefix + "ntracks32" + suffix,
    //                      threeplanes.second.harder_softer_ntracks);
    // treeEvents.push_back(prefix + "nparticles32" + suffix,
    //                      threeplanes.second.harder_softer_nparticles);
    // treeEvents.push_back(prefix + "flavour32" + suffix,
    //                      JetBranch::GetIFNFlavour(threeplanes.second.harder_softer));
    // treeEvents.push_back(prefix + "pTD32" + suffix, JetObservable32.pTD());
    // treeEvents.push_back(prefix + "sigma32" + suffix, width32.sigma);
    // treeEvents.push_back(prefix + "sigma132" + suffix, width32.sigma1);
    // treeEvents.push_back(prefix + "sigma232" + suffix, width32.sigma2);

    // treeEvents.push_back(prefix + "pt41" + suffix, threeplanes.second.softer_harder.pt());
    // treeEvents.push_back(prefix + "eta41" + suffix,
    //                      threeplanes.second.softer_harder.eta());
    // treeEvents.push_back(prefix + "phi41" + suffix,
    //                      threeplanes.second.softer_harder.phi());
    // treeEvents.push_back(prefix + "e41" + suffix, threeplanes.second.softer_harder.e());
    // treeEvents.push_back(prefix + "ntracks41" + suffix,
    //                      threeplanes.second.softer_harder_ntracks);
    // treeEvents.push_back(prefix + "nparticles41" + suffix,
    //                      threeplanes.second.softer_harder_nparticles);
    // treeEvents.push_back(prefix + "flavour41" + suffix,
    //                      JetBranch::GetIFNFlavour(threeplanes.second.softer_harder));
    // treeEvents.push_back(prefix + "pTD41" + suffix, JetObservable41.pTD());
    // treeEvents.push_back(prefix + "sigma41" + suffix, width41.sigma);
    // treeEvents.push_back(prefix + "sigma141" + suffix, width41.sigma1);
    // treeEvents.push_back(prefix + "sigma241" + suffix, width41.sigma2);

    // treeEvents.push_back(prefix + "pt42" + suffix, threeplanes.second.softer_softer.pt());
    // treeEvents.push_back(prefix + "eta42" + suffix,
    //                      threeplanes.second.softer_softer.eta());
    // treeEvents.push_back(prefix + "phi42" + suffix,
    //                      threeplanes.second.softer_softer.phi());
    // treeEvents.push_back(prefix + "e42" + suffix, threeplanes.second.softer_softer.e());
    // treeEvents.push_back(prefix + "ntracks42" + suffix,
    //                      threeplanes.second.softer_softer_ntracks);
    // treeEvents.push_back(prefix + "nparticles42" + suffix,
    //                      threeplanes.second.softer_softer_nparticles);
    // treeEvents.push_back(prefix + "flavour42" + suffix,
    //                      JetBranch::GetIFNFlavour(threeplanes.second.softer_softer));
    // treeEvents.push_back(prefix + "pTD42" + suffix, JetObservable42.pTD());
    // treeEvents.push_back(prefix + "sigma42" + suffix, width42.sigma);
    // treeEvents.push_back(prefix + "sigma142" + suffix, width42.sigma1);
    // treeEvents.push_back(prefix + "sigma242" + suffix, width42.sigma2);

    // treeEvents.push_back(prefix + "pt51" + suffix, threeplanes.third.harder_harder.pt());
    // treeEvents.push_back(prefix + "eta51" + suffix,
    //                      threeplanes.third.harder_harder.eta());
    // treeEvents.push_back(prefix + "phi51" + suffix,
    //                      threeplanes.third.harder_harder.phi());
    // treeEvents.push_back(prefix + "e51" + suffix, threeplanes.third.harder_harder.e());
    // treeEvents.push_back(prefix + "ntracks51" + suffix,
    //                      threeplanes.third.harder_harder_ntracks);
    // treeEvents.push_back(prefix + "nparticles51" + suffix,
    //                      threeplanes.third.harder_harder_nparticles);
    // treeEvents.push_back(prefix + "flavour51" + suffix,
    //                      JetBranch::GetIFNFlavour(threeplanes.third.harder_harder));
    // treeEvents.push_back(prefix + "pTD51" + suffix, JetObservable51.pTD());
    // treeEvents.push_back(prefix + "sigma51" + suffix, width51.sigma);
    // treeEvents.push_back(prefix + "sigma151" + suffix, width51.sigma1);
    // treeEvents.push_back(prefix + "sigma251" + suffix, width51.sigma2);

    // treeEvents.push_back(prefix + "pt52" + suffix, threeplanes.third.harder_softer.pt());
    // treeEvents.push_back(prefix + "eta52" + suffix,
    //                      threeplanes.third.harder_softer.eta());
    // treeEvents.push_back(prefix + "phi52" + suffix,
    //                      threeplanes.third.harder_softer.phi());
    // treeEvents.push_back(prefix + "e52" + suffix, threeplanes.third.harder_softer.e());
    // treeEvents.push_back(prefix + "ntracks52" + suffix,
    //                      threeplanes.third.harder_softer_ntracks);
    // treeEvents.push_back(prefix + "nparticles52" + suffix,
    //                      threeplanes.third.harder_softer_nparticles);
    // treeEvents.push_back(prefix + "flavour52" + suffix,
    //                      JetBranch::GetIFNFlavour(threeplanes.third.harder_softer));
    // treeEvents.push_back(prefix + "pTD52" + suffix, JetObservable52.pTD());
    // treeEvents.push_back(prefix + "sigma52" + suffix, width52.sigma);
    // treeEvents.push_back(prefix + "sigma152" + suffix, width52.sigma1);
    // treeEvents.push_back(prefix + "sigma252" + suffix, width52.sigma2);

    // treeEvents.push_back(prefix + "pt61" + suffix, threeplanes.third.softer_harder.pt());
    // treeEvents.push_back(prefix + "eta61" + suffix,
    //                      threeplanes.third.softer_harder.eta());
    // treeEvents.push_back(prefix + "phi61" + suffix,
    //                      threeplanes.third.softer_harder.phi());
    // treeEvents.push_back(prefix + "e61" + suffix, threeplanes.third.softer_harder.e());
    // treeEvents.push_back(prefix + "ntracks61" + suffix,
    //                      threeplanes.third.softer_harder_ntracks);
    // treeEvents.push_back(prefix + "nparticles61" + suffix,
    //                      threeplanes.third.softer_harder_nparticles);
    // treeEvents.push_back(prefix + "flavour61" + suffix,
    //                      JetBranch::GetIFNFlavour(threeplanes.third.softer_harder));
    // treeEvents.push_back(prefix + "pTD61" + suffix, JetObservable61.pTD());
    // treeEvents.push_back(prefix + "sigma61" + suffix, width61.sigma);
    // treeEvents.push_back(prefix + "sigma161" + suffix, width61.sigma1);
    // treeEvents.push_back(prefix + "sigma261" + suffix, width61.sigma2);

    // treeEvents.push_back(prefix + "pt62" + suffix, threeplanes.third.softer_softer.pt());
    // treeEvents.push_back(prefix + "eta62" + suffix,
    //                      threeplanes.third.softer_softer.eta());
    // treeEvents.push_back(prefix + "phi62" + suffix,
    //                      threeplanes.third.softer_softer.phi());
    // treeEvents.push_back(prefix + "e62" + suffix, threeplanes.third.softer_softer.e());
    // treeEvents.push_back(prefix + "ntracks62" + suffix,
    //                      threeplanes.third.softer_softer_ntracks);
    // treeEvents.push_back(prefix + "nparticles62" + suffix,
    //                      threeplanes.third.softer_softer_nparticles);
    // treeEvents.push_back(prefix + "flavour62" + suffix,
    //                      JetBranch::GetIFNFlavour(threeplanes.third.softer_softer));
    // treeEvents.push_back(prefix + "pTD62" + suffix, JetObservable62.pTD());
    // treeEvents.push_back(prefix + "sigma62" + suffix, width62.sigma);
    // treeEvents.push_back(prefix + "sigma162" + suffix, width62.sigma1);
    // treeEvents.push_back(prefix + "sigma262" + suffix, width62.sigma2);

    treeEvents.push_back(prefix + "init_pdgid1" + suffix,
                         threeplanes.first.harder_init_pdgid);
    treeEvents.push_back(prefix + "init_descri1" + suffix,
                         threeplanes.first.harder_init_descri);
    treeEvents.push_back(prefix + "final_pdgid1" + suffix,
                         threeplanes.first.harder_final_pdgid);
    treeEvents.push_back(prefix + "final_descri1" + suffix,
                         threeplanes.first.harder_final_descri);

    treeEvents.push_back(prefix + "init_pdgid2" + suffix,
                         threeplanes.first.softer_init_pdgid);
    treeEvents.push_back(prefix + "init_descri2" + suffix,
                         threeplanes.first.softer_init_descri);
    treeEvents.push_back(prefix + "final_pdgid2" + suffix,
                         threeplanes.first.softer_final_pdgid);
    treeEvents.push_back(prefix + "final_descri2" + suffix,
                         threeplanes.first.softer_final_descri);

    treeEvents.push_back(prefix + "init_pdgid3" + suffix,
                         threeplanes.second.harder_init_pdgid);
    treeEvents.push_back(prefix + "init_descri3" + suffix,
                         threeplanes.second.harder_init_descri);
    treeEvents.push_back(prefix + "final_pdgid3" + suffix,
                         threeplanes.second.harder_final_pdgid);
    treeEvents.push_back(prefix + "final_descri3" + suffix,
                         threeplanes.second.harder_final_descri);

    treeEvents.push_back(prefix + "init_pdgid4" + suffix,
                         threeplanes.second.softer_init_pdgid);
    treeEvents.push_back(prefix + "init_descri4" + suffix,
                         threeplanes.second.softer_init_descri);
    treeEvents.push_back(prefix + "final_pdgid4" + suffix,
                         threeplanes.second.softer_final_pdgid);
    treeEvents.push_back(prefix + "final_descri4" + suffix,
                         threeplanes.second.softer_final_descri);

    treeEvents.push_back(prefix + "init_pdgid5" + suffix,
                         threeplanes.third.harder_init_pdgid);
    treeEvents.push_back(prefix + "init_descri5" + suffix,
                         threeplanes.third.harder_init_descri);
    treeEvents.push_back(prefix + "final_pdgid5" + suffix,
                         threeplanes.third.harder_final_pdgid);
    treeEvents.push_back(prefix + "final_descri5" + suffix,
                         threeplanes.third.harder_final_descri);

    treeEvents.push_back(prefix + "init_pdgid6" + suffix,
                         threeplanes.third.softer_init_pdgid);
    treeEvents.push_back(prefix + "init_descri6" + suffix,
                         threeplanes.third.softer_init_descri);
    treeEvents.push_back(prefix + "final_pdgid6" + suffix,
                         threeplanes.third.softer_final_pdgid);
    treeEvents.push_back(prefix + "final_descri6" + suffix,
                         threeplanes.third.softer_final_descri);

    // treeEvents.push_back(prefix + "init_pdgid31" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.second.harder_harder, false));
    // treeEvents.push_back(prefix + "init_descri31" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.second.harder_harder, false));
    // treeEvents.push_back(prefix + "final_pdgid31" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.second.harder_harder, true));
    // treeEvents.push_back(prefix + "final_descri31" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.second.harder_harder, true));

    // treeEvents.push_back(prefix + "init_pdgid32" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.second.harder_softer, false));
    // treeEvents.push_back(prefix + "init_descri32" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.second.harder_softer, false));
    // treeEvents.push_back(prefix + "final_pdgid32" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.second.harder_softer, true));
    // treeEvents.push_back(prefix + "final_descri32" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.second.harder_softer, true));

    // treeEvents.push_back(prefix + "init_pdgid41" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.second.softer_harder, false));
    // treeEvents.push_back(prefix + "init_descri41" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.second.softer_harder, false));
    // treeEvents.push_back(prefix + "final_pdgid41" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.second.softer_harder, true));
    // treeEvents.push_back(prefix + "final_descri41" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.second.softer_harder, true));

    // treeEvents.push_back(prefix + "init_pdgid42" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.second.softer_softer, false));
    // treeEvents.push_back(prefix + "init_descri42" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.second.softer_softer, false));
    // treeEvents.push_back(prefix + "final_pdgid42" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.second.softer_softer, true));
    // treeEvents.push_back(prefix + "final_descri42" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.second.softer_softer, true));

    // treeEvents.push_back(prefix + "init_pdgid51" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.third.harder_harder, false));
    // treeEvents.push_back(prefix + "init_descri51" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.third.harder_harder, false));
    // treeEvents.push_back(prefix + "final_pdgid51" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.third.harder_harder, true));
    // treeEvents.push_back(prefix + "final_descri51" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.third.harder_harder, true));

    // treeEvents.push_back(prefix + "init_pdgid52" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.third.harder_softer, false));
    // treeEvents.push_back(prefix + "init_descri52" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.third.harder_softer, false));
    // treeEvents.push_back(prefix + "final_pdgid52" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.third.harder_softer, true));
    // treeEvents.push_back(prefix + "final_descri52" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.third.harder_softer, true));

    // treeEvents.push_back(prefix + "init_pdgid61" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.third.softer_harder, false));
    // treeEvents.push_back(prefix + "init_descri61" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.third.softer_harder, false));
    // treeEvents.push_back(prefix + "final_pdgid61" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.third.softer_harder, true));
    // treeEvents.push_back(prefix + "final_descri61" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.third.softer_harder, true));

    // treeEvents.push_back(prefix + "init_pdgid62" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.third.softer_softer, false));
    // treeEvents.push_back(prefix + "init_descri62" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.third.softer_softer, false));
    // treeEvents.push_back(prefix + "final_pdgid62" + suffix,
    //                      JetBranch::GetIFNPdgid(threeplanes.third.softer_softer, true));
    // treeEvents.push_back(prefix + "final_descri62" + suffix,
    //                      JetBranch::GetIFNDescri(threeplanes.third.softer_softer, true));

    treeEvents.push_back(prefix + "nbHadrons1" + suffix,
                         threeplanes.first.harder_ghostb);
    treeEvents.push_back(prefix + "nbHadrons2" + suffix,
                         threeplanes.first.softer_ghostb);
    treeEvents.push_back(prefix + "nbHadrons3" + suffix,
                         threeplanes.second.harder_ghostb);
    treeEvents.push_back(prefix + "nbHadrons4" + suffix,
                         threeplanes.second.softer_ghostb);
    treeEvents.push_back(prefix + "nbHadrons5" + suffix,
                         threeplanes.third.harder_ghostb);
    treeEvents.push_back(prefix + "nbHadrons6" + suffix,
                         threeplanes.third.softer_ghostb);

    treeEvents.push_back(prefix + "ncHadrons1" + suffix,
                         threeplanes.first.harder_ghostc);
    treeEvents.push_back(prefix + "ncHadrons2" + suffix,
                         threeplanes.first.softer_ghostc);
    treeEvents.push_back(prefix + "ncHadrons3" + suffix,
                         threeplanes.second.harder_ghostc);
    treeEvents.push_back(prefix + "ncHadrons4" + suffix,
                         threeplanes.second.softer_ghostc);
    treeEvents.push_back(prefix + "ncHadrons5" + suffix,
                         threeplanes.third.harder_ghostc);
    treeEvents.push_back(prefix + "ncHadrons6" + suffix,
                         threeplanes.third.softer_ghostc);

    treeEvents.push_back(prefix + "z1" + suffix, threeplanes.first.z);
    treeEvents.push_back(prefix + "kt1" + suffix, threeplanes.first.kt);
    treeEvents.push_back(prefix + "theta1" + suffix, threeplanes.first.theta);
    treeEvents.push_back(prefix + "deltaR1" + suffix, threeplanes.first.deltaR);
    treeEvents.push_back(prefix + "eec1" + suffix,
                         (threeplanes.first.harder.e() * threeplanes.first.softer.e()) * 1.0 /
                             threeplanes.Q / threeplanes.Q);
    treeEvents.push_back(prefix + "type1" + suffix, JetObservable::determineType(threeplanes.first.softer_flav, threeplanes.first.harder_flav));

    treeEvents.push_back(prefix + "z2" + suffix, threeplanes.second.z);
    treeEvents.push_back(prefix + "kt2" + suffix, threeplanes.second.kt);
    treeEvents.push_back(prefix + "theta2" + suffix, threeplanes.second.theta);
    treeEvents.push_back(prefix + "deltaR2" + suffix, threeplanes.second.deltaR);
    treeEvents.push_back(
        prefix + "eec2" + suffix,
        (threeplanes.second.harder.e() * threeplanes.second.softer.e()) * 1.0 /
            threeplanes.Q / threeplanes.Q);
    treeEvents.push_back(prefix + "type2" + suffix, JetObservable::determineType(threeplanes.second.softer_flav, threeplanes.second.harder_flav));

    treeEvents.push_back(prefix + "z3" + suffix, threeplanes.third.z);
    treeEvents.push_back(prefix + "kt3" + suffix, threeplanes.third.kt);
    treeEvents.push_back(prefix + "theta3" + suffix, threeplanes.third.theta);
    treeEvents.push_back(prefix + "deltaR3" + suffix, threeplanes.third.deltaR);
    treeEvents.push_back(
        prefix + "eec3" + suffix,
        (threeplanes.third.harder.e() * threeplanes.third.softer.e()) * 1.0 /
            threeplanes.Q / threeplanes.Q);
    treeEvents.push_back(prefix + "type3" + suffix, JetObservable::determineType(threeplanes.third.softer_flav, threeplanes.third.harder_flav));

    treeEvents.push_back(prefix + "n1x" + suffix, threeplanes.first.n.x());
    treeEvents.push_back(prefix + "n1y" + suffix, threeplanes.first.n.y());
    treeEvents.push_back(prefix + "n1z" + suffix, threeplanes.first.n.z());
    treeEvents.push_back(prefix + "n2x" + suffix, threeplanes.second.n.x());
    treeEvents.push_back(prefix + "n2y" + suffix, threeplanes.second.n.y());
    treeEvents.push_back(prefix + "n2z" + suffix, threeplanes.second.n.z());
    treeEvents.push_back(prefix + "n3x" + suffix, threeplanes.third.n.x());
    treeEvents.push_back(prefix + "n3y" + suffix, threeplanes.third.n.y());
    treeEvents.push_back(prefix + "n3z" + suffix, threeplanes.third.n.z());

    treeEvents.push_back(prefix + "dPhi12_X2" + suffix, planetheta2.dphi12_X);
    treeEvents.push_back(prefix + "Theta2" + suffix, planetheta2.theta);
    treeEvents.push_back(prefix + "Theta22" + suffix, planetheta2.theta2);
    treeEvents.push_back(prefix + "dPhi12_X3" + suffix, planetheta3.dphi12_X);
    treeEvents.push_back(prefix + "Theta3" + suffix, planetheta3.theta);
    treeEvents.push_back(prefix + "Theta23" + suffix, planetheta3.theta2);
    treeEvents.push_back(prefix + "Phi2" + suffix, threeplanes.deltaPhi);
    treeEvents.push_back(prefix + "Phi3" + suffix, threeplanes.deltaPhi2);

    treeEvents.push_back(prefix + "jpt" + suffix, threeplanes.first.initJet.Pt());
    treeEvents.push_back(prefix + "jeta" + suffix,
                         threeplanes.first.initJet.Eta());
    treeEvents.push_back(prefix + "jphi" + suffix,
                         threeplanes.first.initJet.Phi());
    treeEvents.push_back(prefix + "je" + suffix, threeplanes.first.initJet.E());

    if (SampleType.find("CMS") != std::string::npos)
    {
      treeEvents.push_back(prefix + "subjetbpt1" + suffix, threeplanes.bchadrons.bhadron1.pt());
      treeEvents.push_back(prefix + "subjetbeta1" + suffix,
                           threeplanes.bchadrons.bhadron1.eta());
      treeEvents.push_back(prefix + "subjetbphi1" + suffix,
                           threeplanes.bchadrons.bhadron1.phi());
      treeEvents.push_back(prefix + "subjetbe1" + suffix, threeplanes.bchadrons.bhadron1.e());

      treeEvents.push_back(prefix + "subjetbpt2" + suffix, threeplanes.bchadrons.bhadron2.pt());
      treeEvents.push_back(prefix + "subjetbeta2" + suffix,
                           threeplanes.bchadrons.bhadron2.eta());
      treeEvents.push_back(prefix + "subjetbphi2" + suffix,
                           threeplanes.bchadrons.bhadron2.phi());
      treeEvents.push_back(prefix + "subjetbe2" + suffix, threeplanes.bchadrons.bhadron2.e());

      treeEvents.push_back(prefix + "subjetbpt3" + suffix, threeplanes.bchadrons.bhadron3.pt());
      treeEvents.push_back(prefix + "subjetbeta3" + suffix,
                           threeplanes.bchadrons.bhadron3.eta());
      treeEvents.push_back(prefix + "subjetbphi3" + suffix,
                           threeplanes.bchadrons.bhadron3.phi());
      treeEvents.push_back(prefix + "subjetbe3" + suffix, threeplanes.bchadrons.bhadron3.e());

      treeEvents.push_back(prefix + "subjetbpt4" + suffix, threeplanes.bchadrons.bhadron4.pt());
      treeEvents.push_back(prefix + "subjetbeta4" + suffix,
                           threeplanes.bchadrons.bhadron4.eta());
      treeEvents.push_back(prefix + "subjetbphi4" + suffix,
                           threeplanes.bchadrons.bhadron4.phi());
      treeEvents.push_back(prefix + "subjetbe4" + suffix, threeplanes.bchadrons.bhadron4.e());

      treeEvents.push_back(prefix + "subjetbPhi" + suffix, threeplanes.bchadrons.bhadronphi);
      treeEvents.push_back(prefix + "subjetbdPhi12_X" + suffix,
                           threeplanes.bchadrons.bhadronboostphi);

      treeEvents.push_back(prefix + "subjetcpt1" + suffix, threeplanes.bchadrons.chadron1.pt());
      treeEvents.push_back(prefix + "subjetceta1" + suffix,
                           threeplanes.bchadrons.chadron1.eta());
      treeEvents.push_back(prefix + "subjetcphi1" + suffix,
                           threeplanes.bchadrons.chadron1.phi());
      treeEvents.push_back(prefix + "subjetce1" + suffix, threeplanes.bchadrons.chadron1.e());

      treeEvents.push_back(prefix + "subjetcpt2" + suffix, threeplanes.bchadrons.chadron2.pt());
      treeEvents.push_back(prefix + "subjetceta2" + suffix,
                           threeplanes.bchadrons.chadron2.eta());
      treeEvents.push_back(prefix + "subjetcphi2" + suffix,
                           threeplanes.bchadrons.chadron2.phi());
      treeEvents.push_back(prefix + "subjetce2" + suffix, threeplanes.bchadrons.chadron2.e());

      treeEvents.push_back(prefix + "subjetcpt3" + suffix, threeplanes.bchadrons.chadron3.pt());
      treeEvents.push_back(prefix + "subjetceta3" + suffix,
                           threeplanes.bchadrons.chadron3.eta());
      treeEvents.push_back(prefix + "subjetcphi3" + suffix,
                           threeplanes.bchadrons.chadron3.phi());
      treeEvents.push_back(prefix + "subjetce3" + suffix, threeplanes.bchadrons.chadron3.e());

      treeEvents.push_back(prefix + "subjetcpt4" + suffix, threeplanes.bchadrons.chadron4.pt());
      treeEvents.push_back(prefix + "subjetceta4" + suffix,
                           threeplanes.bchadrons.chadron4.eta());
      treeEvents.push_back(prefix + "subjetcphi4" + suffix,
                           threeplanes.bchadrons.chadron4.phi());
      treeEvents.push_back(prefix + "subjetce4" + suffix, threeplanes.bchadrons.chadron4.e());

      treeEvents.push_back(prefix + "subjetcPhi" + suffix, threeplanes.bchadrons.chadronphi);
      treeEvents.push_back(prefix + "subjetcdPhi12_X" + suffix,
                           threeplanes.bchadrons.chadronboostphi);

      treeEvents.push_back(prefix + "nbhadrons" + suffix, threeplanes.bchadrons.nbhadrons);
      treeEvents.push_back(prefix + "nchadrons" + suffix,
                           threeplanes.bchadrons.nchadrons);
    }
    if (SaveParticles)
    {
      for (int ii = 0; ii < threeplanes.first.harder_constituents_info.size();
           ii++)
      {
        auto particle = threeplanes.first.harder_constituents_info.at(ii);
        treeEvents.push_back(prefix + "particle1_pt" + suffix, particle.pt);
        treeEvents.push_back(prefix + "particle1_eta" + suffix, particle.eta);
        treeEvents.push_back(prefix + "particle1_phi" + suffix, particle.phi);
        treeEvents.push_back(prefix + "particle1_e" + suffix, particle.e);
        treeEvents.push_back(prefix + "particle1_charge" + suffix,
                             particle.charge);
        treeEvents.push_back(prefix + "particle1_pid" + suffix, particle.pdgid);
        treeEvents.push_back(prefix + "particle1_jetid" + suffix, jetindex);
      }
      for (int ii = 0; ii < threeplanes.first.softer_constituents_info.size();
           ii++)
      {
        auto particle = threeplanes.first.softer_constituents_info.at(ii);
        treeEvents.push_back(prefix + "particle2_pt" + suffix, particle.pt);
        treeEvents.push_back(prefix + "particle2_eta" + suffix, particle.eta);
        treeEvents.push_back(prefix + "particle2_phi" + suffix, particle.phi);
        treeEvents.push_back(prefix + "particle2_e" + suffix, particle.e);
        treeEvents.push_back(prefix + "particle2_charge" + suffix,
                             particle.charge);
        treeEvents.push_back(prefix + "particle2_pid" + suffix, particle.pdgid);
        treeEvents.push_back(prefix + "particle2_jetid" + suffix, jetindex);
      }

      for (int ii = 0; ii < threeplanes.second.harder_constituents_info.size();
           ii++)
      {
        auto particle = threeplanes.second.harder_constituents_info.at(ii);
        treeEvents.push_back(prefix + "particle3_pt" + suffix, particle.pt);
        treeEvents.push_back(prefix + "particle3_eta" + suffix, particle.eta);
        treeEvents.push_back(prefix + "particle3_phi" + suffix, particle.phi);
        treeEvents.push_back(prefix + "particle3_e" + suffix, particle.e);
        treeEvents.push_back(prefix + "particle3_charge" + suffix,
                             particle.charge);
        treeEvents.push_back(prefix + "particle3_pid" + suffix, particle.pdgid);
        treeEvents.push_back(prefix + "particle3_jetid" + suffix, jetindex);
      }
      for (int ii = 0; ii < threeplanes.second.softer_constituents_info.size();
           ii++)
      {
        auto particle = threeplanes.second.softer_constituents_info.at(ii);
        treeEvents.push_back(prefix + "particle4_pt" + suffix, particle.pt);
        treeEvents.push_back(prefix + "particle4_eta" + suffix, particle.eta);
        treeEvents.push_back(prefix + "particle4_phi" + suffix, particle.phi);
        treeEvents.push_back(prefix + "particle4_e" + suffix, particle.e);
        treeEvents.push_back(prefix + "particle4_charge" + suffix,
                             particle.charge);
        treeEvents.push_back(prefix + "particle4_pid" + suffix, particle.pdgid);
        treeEvents.push_back(prefix + "particle4_jetid" + suffix, jetindex);
      }

      for (int ii = 0; ii < threeplanes.third.harder_constituents_info.size();
           ii++)
      {
        auto particle = threeplanes.third.harder_constituents_info.at(ii);
        treeEvents.push_back(prefix + "particle5_pt" + suffix, particle.pt);
        treeEvents.push_back(prefix + "particle5_eta" + suffix, particle.eta);
        treeEvents.push_back(prefix + "particle5_phi" + suffix, particle.phi);
        treeEvents.push_back(prefix + "particle5_e" + suffix, particle.e);
        treeEvents.push_back(prefix + "particle5_charge" + suffix,
                             particle.charge);
        treeEvents.push_back(prefix + "particle5_pid" + suffix, particle.pdgid);
        treeEvents.push_back(prefix + "particle5_jetid" + suffix, jetindex);
      }
      for (int ii = 0; ii < threeplanes.third.softer_constituents_info.size();
           ii++)
      {
        auto particle = threeplanes.third.softer_constituents_info.at(ii);
        treeEvents.push_back(prefix + "particle6_pt" + suffix, particle.pt);
        treeEvents.push_back(prefix + "particle6_eta" + suffix, particle.eta);
        treeEvents.push_back(prefix + "particle6_phi" + suffix, particle.phi);
        treeEvents.push_back(prefix + "particle6_e" + suffix, particle.e);
        treeEvents.push_back(prefix + "particle6_charge" + suffix,
                             particle.charge);
        treeEvents.push_back(prefix + "particle6_pid" + suffix, particle.pdgid);
        treeEvents.push_back(prefix + "particle6_jetid" + suffix, jetindex);
      }
    }
  }
};
