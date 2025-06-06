#pragma once
#include "TLorentzVector.h"
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <random>
#include <tuple>
#include <unordered_map>
#include <vector>
#include "Matching.cc"
#include "ParticleInfo.h"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include "fastjet/contrib/IFNPlugin.hh"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "SpinObservable.h"

class Matching;
using namespace fastjet;
using namespace fastjet::contrib;
using namespace std;
class JetBranch
{
public:
  struct BChadronPlanes
  {
    PseudoJet bhadron1;
    PseudoJet bhadron2;
    PseudoJet bhadron3;
    PseudoJet bhadron4;
    double bhadronphi;
    double bhadronboostphi;
    PseudoJet chadron1;
    PseudoJet chadron2;
    PseudoJet chadron3;
    PseudoJet chadron4;
    double chadronphi;
    double chadronboostphi;
    int nbhadrons;
    int nchadrons;
  };
  struct PlaneVariables
  {
    PseudoJet harder;
    PseudoJet softer;
    PseudoJet harder_harder;
    PseudoJet harder_softer;
    PseudoJet softer_harder;
    PseudoJet softer_softer;
    TLorentzVector harder_vect;
    TLorentzVector softer_vect;
    double theta;
    double deltaR;
    double z;
    double kt;
    TVector3 n;
    vector<PseudoJet> harder_constituents;
    vector<PseudoJet> softer_constituents;
    vector<PseudoJet> harder_harder_constituents;
    vector<PseudoJet> harder_softer_constituents;
    vector<PseudoJet> softer_harder_constituents;
    vector<PseudoJet> softer_softer_constituents;
    int harder_ntracks;
    int softer_ntracks;
    int harder_nparticles;
    int softer_nparticles;
    int harder_harder_ntracks;
    int softer_harder_ntracks;
    int harder_harder_nparticles;
    int softer_harder_nparticles;
    int harder_softer_ntracks;
    int softer_softer_ntracks;
    int harder_softer_nparticles;
    int softer_softer_nparticles;
    vector<ParticleInfo> harder_constituents_info;
    vector<ParticleInfo> softer_constituents_info;
    vector<ParticleInfo> harder_harder_constituents_info;
    vector<ParticleInfo> harder_softer_constituents_info;
    vector<ParticleInfo> softer_harder_constituents_info;
    vector<ParticleInfo> softer_softer_constituents_info;
    int softer_flav;
    int harder_flav;
    int softer_init_pdgid;
    int softer_init_descri;
    int softer_final_pdgid;
    int softer_final_descri;
    int harder_init_pdgid;
    int harder_init_descri;
    int harder_final_pdgid;
    int harder_final_descri;
    int harder_ghostb = 0;
    int harder_ghostc = 0;
    int softer_ghostb = 0;
    int softer_ghostc = 0;
    void GetVector(vector<ParticleInfo> particlesinfo)
    {
      int harder_ghostb_temp = 0;
      int harder_ghostc_temp = 0;
      int softer_ghostb_temp = 0;
      int softer_ghostc_temp = 0;
      harder_vect = JetBranch::PseudoJetToTLorentzVector(harder);
      softer_vect = JetBranch::PseudoJetToTLorentzVector(softer);
      GetPlaneVector();
      harder_ntracks = GetNtracks(harder_constituents, particlesinfo);
      softer_ntracks = GetNtracks(softer_constituents, particlesinfo);
      harder_harder_ntracks = GetNtracks(harder_harder_constituents, particlesinfo);
      softer_harder_ntracks = GetNtracks(softer_harder_constituents, particlesinfo);
      harder_softer_ntracks = GetNtracks(harder_softer_constituents, particlesinfo);
      softer_softer_ntracks = GetNtracks(softer_softer_constituents, particlesinfo);
      for (int i = 0; i < harder_constituents.size(); i++)
      {
        auto particle = harder_constituents.at(i);
        harder_constituents_info.push_back(ParticleInfo(
            particlesinfo.at(particle.user_index()).pdgid,
            particlesinfo.at(particle.user_index()).charge, particle.pt(),
            particle.eta(), particle.phi(), particle.e()));
        if ((abs(particlesinfo.at(particle.user_index()).pdgid / 100 % 10) == 5 ||
             abs(particlesinfo.at(particle.user_index()).pdgid / 1000 % 10) == 5))
          harder_ghostb_temp++;
        else if ((abs(particlesinfo.at(particle.user_index()).pdgid / 100 % 10) == 4 ||
                  abs(particlesinfo.at(particle.user_index()).pdgid / 1000 % 10) == 4))
          harder_ghostc_temp++;
      }
      for (int i = 0; i < softer_constituents.size(); i++)
      {
        auto particle = softer_constituents.at(i);
        softer_constituents_info.push_back(ParticleInfo(
            particlesinfo.at(particle.user_index()).pdgid,
            particlesinfo.at(particle.user_index()).charge, particle.pt(),
            particle.eta(), particle.phi(), particle.e()));
        if ((abs(particlesinfo.at(particle.user_index()).pdgid / 100 % 10) == 5 ||
             abs(particlesinfo.at(particle.user_index()).pdgid / 1000 % 10) == 5))
          softer_ghostb_temp++;
        else if ((abs(particlesinfo.at(particle.user_index()).pdgid / 100 % 10) == 4 ||
                  abs(particlesinfo.at(particle.user_index()).pdgid / 1000 % 10) == 4))
          softer_ghostc_temp++;
      }
      harder_ghostb = harder_ghostb_temp;
      harder_ghostc = harder_ghostc_temp;
      softer_ghostb = softer_ghostb_temp;
      softer_ghostc = softer_ghostc_temp;
      harder_nparticles = harder_constituents.size() - harder_ghostb - harder_ghostc;
      softer_nparticles = softer_constituents.size() - softer_ghostb - softer_ghostc;
      harder_harder_nparticles = harder_harder_constituents.size();
      softer_harder_nparticles = softer_harder_constituents.size();
      harder_softer_nparticles = harder_softer_constituents.size();
      softer_softer_nparticles = softer_softer_constituents.size();
      for (int i = 0; i < harder_harder_constituents.size(); i++)
      {
        auto particle = harder_harder_constituents.at(i);
        harder_harder_constituents_info.push_back(ParticleInfo(
            particlesinfo.at(particle.user_index()).pdgid,
            particlesinfo.at(particle.user_index()).charge, particle.pt(),
            particle.eta(), particle.phi(), particle.e()));
      }
      for (int i = 0; i < harder_softer_constituents.size(); i++)
      {
        auto particle = harder_softer_constituents.at(i);
        harder_softer_constituents_info.push_back(ParticleInfo(
            particlesinfo.at(particle.user_index()).pdgid,
            particlesinfo.at(particle.user_index()).charge, particle.pt(),
            particle.eta(), particle.phi(), particle.e()));
      }
      for (int i = 0; i < softer_harder_constituents.size(); i++)
      {
        auto particle = softer_harder_constituents.at(i);
        softer_harder_constituents_info.push_back(ParticleInfo(
            particlesinfo.at(particle.user_index()).pdgid,
            particlesinfo.at(particle.user_index()).charge, particle.pt(),
            particle.eta(), particle.phi(), particle.e()));
      }
      for (int i = 0; i < softer_softer_constituents.size(); i++)
      {
        auto particle = softer_softer_constituents.at(i);
        softer_softer_constituents_info.push_back(ParticleInfo(
            particlesinfo.at(particle.user_index()).pdgid,
            particlesinfo.at(particle.user_index()).charge, particle.pt(),
            particle.eta(), particle.phi(), particle.e()));
      }
    }
    void GetPlaneVector()
    {
      n = harder_vect.Vect().Cross(softer_vect.Vect());
    }

    bool isqq = false;
    bool isgg = false;
    TLorentzVector initJet;
  };
  struct twoplanes
  {
    JetBranch::PlaneVariables first;
    JetBranch::PlaneVariables second;
    double deltaPhi;
    double Q = 0;
  };
  struct pseudothreeplanes
  {
    JetBranch::PlaneVariables first;
    JetBranch::PlaneVariables second;
    JetBranch::PlaneVariables third;
    JetBranch::BChadronPlanes bchadrons;
    double deltaPhi = -10;
    double deltaPhi2 = -10;
    double Q = 0;
    bool ismatched = false;
    int firstindex = 0;
    pseudothreeplanes() = default;
    pseudothreeplanes(const JetBranch::PlaneVariables &f, const JetBranch::PlaneVariables &s, const JetBranch::PlaneVariables &t,
                      double dPhi, double dPhi2, double q, bool _ismatched, int _firstindex, JetBranch::BChadronPlanes _bchadrons)
        : first(f), second(s), third(t), deltaPhi(dPhi), deltaPhi2(dPhi2), Q(q), ismatched(_ismatched), firstindex(_firstindex), bchadrons(_bchadrons) {}
  };

  struct threeplanes
  {
    JetBranch::PlaneVariables first;
    JetBranch::PlaneVariables second;
    JetBranch::PlaneVariables third;
    JetBranch::BChadronPlanes bchadrons;
    double deltaPhi;
    double deltaPhi2;
    double Q = 0;
    pseudothreeplanes matchedthreeplanes;
    bool ismatched = false;
    int firstindex = 0;
    threeplanes() = default;
    pseudothreeplanes Getpseudothreeplanes() const
    {
      double t_deltaPhi = deltaPhi;
      double t_deltaPhi2 = deltaPhi2;
      if (second.softer_nparticles == 0 || second.harder_nparticles == 0)
        t_deltaPhi = -10;
      if (third.softer_nparticles == 0 || third.harder_nparticles == 0)
        t_deltaPhi2 = -10;
      return pseudothreeplanes(first, second, third, t_deltaPhi, t_deltaPhi2, Q, ismatched, firstindex, bchadrons);
    }
  };
  static double DeltaPhi(PlaneVariables plane1, PlaneVariables plane2)
  {
    if (plane1.harder_nparticles == 0 || plane1.softer_nparticles == 0 ||
        plane2.harder_nparticles == 0 || plane2.softer_nparticles == 0)
      return -10;
    TVector3 n1 = plane1.n;
    TVector3 n2 = plane2.n;
    n1 = n1.Unit();
    n2 = n2.Unit();
    double cos_angle = abs(n1.Dot(n2));
    if ((n1.Cross(n2)).Dot(plane2.harder_vect.Vect()) < 0)
    {
      cos_angle = -1.0 * cos_angle;
    }
    double angle_radians = TMath::ACos(cos_angle);
    return angle_radians;
  }

  static TLorentzVector PseudoJetToTLorentzVector(fastjet::PseudoJet &pj)
  {
    return TLorentzVector(pj.px(), pj.py(), pj.pz(), pj.e());
  }

  static twoplanes findPrimaryAndSecondaryJets(PseudoJet &j0, double z1cut, double z2cut, double kt1cut, double kt2cut, std::vector<ParticleInfo> particlesinfo, bool issecondsoft = true)
  {
    PseudoJet j1, j2, j3, j4, jinit;
    PlaneVariables maxPrimary = {
        PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), TLorentzVector(), TLorentzVector(), 0.0, 0.0,
        0.0, 0.0, TVector3()};
    PlaneVariables maxSecondary = {
        PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), TLorentzVector(), TLorentzVector(), 0.0, 0.0,
        0.0, 0.0, TVector3()};
    jinit = j0;

    while (j0.has_parents(j1, j2))
    {

      PseudoJet primaryHarder = (j1.modp() > j2.modp()) ? j1 : j2;
      PseudoJet primarySofter = (j1.modp() > j2.modp()) ? j2 : j1;
      PseudoJet primaryHarder_harder, primaryHarder_softer, j11, j12;
      if (primaryHarder.has_parents(j11, j12))
      {
        primaryHarder_harder = (j11.modp() > j12.modp()) ? j11 : j12;
        primaryHarder_softer = (j11.modp() > j12.modp()) ? j12 : j11;
      }
      PseudoJet primarySofter_harder, primarySofter_softer, j21, j22;
      if (primarySofter.has_parents(j21, j22))
      {
        primarySofter_harder = (j21.modp() > j22.modp()) ? j21 : j22;
        primarySofter_softer = (j21.modp() > j22.modp()) ? j22 : j21;
      }
      double primaryTheta = primaryHarder.delta_phi_to(primarySofter);
      double primaryDeltaR = primaryHarder.delta_R(primarySofter);
      double primaryZ =
          primarySofter.modp() / (primarySofter.modp() + primaryHarder.modp());
      double primaryKt = primarySofter.modp() * primaryDeltaR;
      // double primaryKt = primarySofter.modp() * sin(primaryTheta);
      // Update maxPrimary if this kt is larger
      if (primaryKt > maxPrimary.kt && primaryZ > z1cut && primaryKt < kt1cut)
      {
        maxPrimary = {primaryHarder,
                      primarySofter,
                      primaryHarder_harder,
                      primaryHarder_softer,
                      primarySofter_harder,
                      primarySofter_softer,
                      PseudoJetToTLorentzVector(primaryHarder),
                      PseudoJetToTLorentzVector(primarySofter),
                      primaryTheta,
                      primaryDeltaR,
                      primaryZ,
                      primaryKt,
                      TVector3(),
                      primaryHarder.constituents(),
                      primarySofter.constituents(),
                      primaryHarder_harder.constituents(),
                      primaryHarder_softer.constituents(),
                      primarySofter_harder.constituents(),
                      primarySofter_softer.constituents()};
        maxSecondary = {
            PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), TLorentzVector(), TLorentzVector(), 0.0, 0.0,
            0.0, 0.0, TVector3()};
        PseudoJet secondplanejet = primaryHarder;
        if (issecondsoft)
          secondplanejet = primarySofter;
        while (secondplanejet.has_parents(j3, j4))
        {
          PseudoJet secondaryHarder = (j3.modp() > j4.modp()) ? j3 : j4;
          PseudoJet secondarySofter = (j3.modp() > j4.modp()) ? j4 : j3;
          PseudoJet secondaryHarder_harder, secondaryHarder_softer, j31, j32;
          if (secondaryHarder.has_parents(j31, j32))
          {
            secondaryHarder_harder = (j31.modp() > j32.modp()) ? j31 : j32;
            secondaryHarder_softer = (j31.modp() > j32.modp()) ? j32 : j31;
          }
          PseudoJet secondarySofter_harder, secondarySofter_softer, j41, j42;
          if (secondarySofter.has_parents(j41, j42))
          {
            secondarySofter_harder = (j41.modp() > j42.modp()) ? j41 : j42;
            secondarySofter_softer = (j41.modp() > j42.modp()) ? j42 : j41;
          }
          double secondaryTheta = secondaryHarder.delta_phi_to(secondarySofter);
          double secondaryDeltaR =
              secondaryHarder.delta_R(secondarySofter);
          double secondaryZ = secondarySofter.modp() /
                              (secondarySofter.modp() + secondaryHarder.modp());
          double secondaryKt = secondarySofter.modp() * secondaryDeltaR;
          if (secondaryKt > maxSecondary.kt && secondaryZ > z2cut && secondaryKt > kt2cut)
          {
            maxSecondary = {secondaryHarder,
                            secondarySofter,
                            secondaryHarder_harder,
                            secondaryHarder_softer,
                            secondarySofter_harder,
                            secondarySofter_softer,
                            PseudoJetToTLorentzVector(secondaryHarder),
                            PseudoJetToTLorentzVector(secondarySofter),
                            secondaryTheta,
                            secondaryDeltaR,
                            secondaryZ,
                            secondaryKt,
                            TVector3(),
                            secondaryHarder.constituents(),
                            secondarySofter.constituents(),
                            secondaryHarder_harder.constituents(),
                            secondaryHarder_softer.constituents(),
                            secondarySofter_harder.constituents(),
                            secondarySofter_softer.constituents()};
          }
          secondplanejet = secondaryHarder;
        }
      }

      j0 = primaryHarder; // Move to the next iteration with the harder subjet
    }

    /// Add plane vector
    maxPrimary.GetVector(particlesinfo);
    maxSecondary.GetVector(particlesinfo);
    twoplanes twoplanes_out;
    maxPrimary.initJet.SetPtEtaPhiE(jinit.pt(), jinit.eta(), jinit.phi(), jinit.e());
    maxSecondary.initJet.SetPtEtaPhiE(jinit.pt(), jinit.eta(), jinit.phi(), jinit.e());
    twoplanes_out.first = maxPrimary;
    twoplanes_out.second = maxSecondary;
    twoplanes_out.deltaPhi = JetBranch::DeltaPhi(maxPrimary, maxSecondary);
    return twoplanes_out;
  }
  static BChadronPlanes findBhadronPlanes(PseudoJet &j0, std::vector<ParticleInfo> particlesinfo)
  {
    PseudoJet j1, j2, j3, j4, jinit;
    BChadronPlanes planes = {PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), -10, -10,
                             PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), -10, -10,
                             0, 0};
    jinit = j0;
    int bnum = GetGhostNumber(jinit, particlesinfo, 5);
    int cnum = GetGhostNumber(jinit, particlesinfo, 4);
    int index = 0;
    planes.nbhadrons = bnum;
    planes.nchadrons = cnum;
    if (bnum == 2)
    {
      while (jinit.has_parents(j1, j2))
      {
        index++;
        int bnum_j1 = GetGhostNumber(j1, particlesinfo, 5);
        int bnum_j2 = GetGhostNumber(j2, particlesinfo, 5);
        // std::cout << bnum_j1 << " " << bnum_j2 << std::endl;
        if ((bnum_j1 == 1 && bnum_j2 == 1))
        {
          PseudoJet j_partner;
          if (jinit.has_partner(j_partner))
          {
            SpinObservable spin(PseudoJetToTLorentzVector(jinit), PseudoJetToTLorentzVector(j_partner),
                                PseudoJetToTLorentzVector(j1), PseudoJetToTLorentzVector(j2));
            planes.bhadron1 = jinit;
            planes.bhadron2 = j_partner;
            planes.bhadron3 = j1;
            planes.bhadron4 = j2;
            planes.bhadronphi = spin.GetPhi();
            planes.bhadronboostphi = spin.GetPlaneTheta().dphi12_X;
          }
          else
          {
            planes.bhadronphi = -5;
            planes.bhadronboostphi = -5;
          }
          break;
        }
        else if (bnum_j1 == 2 && bnum_j2 == 0)
        {
          jinit = j1;
        }
        else if (bnum_j1 == 0 && bnum_j2 == 2)
        {
          jinit = j2;
        }
        else
        {
          std::cout << "Error: bnum not right" << std::endl;
        }
      }
    }

    jinit = j0;
    if (cnum == 2)
    {
      while (jinit.has_parents(j1, j2))
      {
        index++;
        int cnum_j1 = GetGhostNumber(j1, particlesinfo, 4);
        int cnum_j2 = GetGhostNumber(j2, particlesinfo, 4);
        if (cnum_j1 == 1 && cnum_j2 == 1)
        {
          PseudoJet j_partner;
          if (jinit.has_partner(j_partner))
          {
            SpinObservable spin(PseudoJetToTLorentzVector(jinit), PseudoJetToTLorentzVector(j_partner),
                                PseudoJetToTLorentzVector(j1), PseudoJetToTLorentzVector(j2));
            planes.chadron1 = jinit;
            planes.chadron2 = j_partner;
            planes.chadron3 = j1;
            planes.chadron4 = j2;
            planes.chadronphi = spin.GetPhi();
            planes.chadronboostphi = spin.GetPlaneTheta().dphi12_X;
          }
          else
          {
            planes.chadronphi = -5;
            planes.chadronboostphi = -5;
          }
          break;
        }
        else if (cnum_j1 == 2 && cnum_j2 == 0)
        {
          jinit = j1;
        }
        else if (cnum_j1 == 0 && cnum_j2 == 2)
        {
          jinit = j2;
        }
        else
        {
          std::cout << "Error: cnum not right" << std::endl;
        }
      }
    }

    return planes;
  }
  static threeplanes findPrimarySecondaryAndThirdaryJets(PseudoJet &j0, double z1cut, double z2cut, double kt1cut, double kt2cut, std::vector<ParticleInfo> particlesinfo, bool issecondsoft = true)
  {
    PseudoJet j1, j2, j3, j4, jinit;
    PlaneVariables maxPrimary = {
        PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), TLorentzVector(), TLorentzVector(), 0.0, 0.0,
        0.0, 0.0, TVector3()};
    PlaneVariables maxSecondary = {
        PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), TLorentzVector(), TLorentzVector(), 0.0, 0.0,
        0.0, 0.0, TVector3()};
    PlaneVariables maxThirdary = {
        PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), TLorentzVector(), TLorentzVector(), 0.0, 0.0,
        0.0, 0.0, TVector3()};
    jinit = j0;
    std::vector<bool> isselects;
    while (j0.has_parents(j1, j2))
    {
      PseudoJet primaryHarder = (j1.modp() > j2.modp()) ? j1 : j2;
      PseudoJet primarySofter = (j1.modp() > j2.modp()) ? j2 : j1;
      PseudoJet primaryHarder_harder, primaryHarder_softer, j11, j12;
      if (primaryHarder.has_parents(j11, j12))
      {
        primaryHarder_harder = (j11.modp() > j12.modp()) ? j11 : j12;
        primaryHarder_softer = (j11.modp() > j12.modp()) ? j12 : j11;
      }
      PseudoJet primarySofter_harder, primarySofter_softer, j21, j22;
      if (primarySofter.has_parents(j21, j22))
      {
        primarySofter_harder = (j21.modp() > j22.modp()) ? j21 : j22;
        primarySofter_softer = (j21.modp() > j22.modp()) ? j22 : j21;
      }
      double primaryTheta = primaryHarder.delta_phi_to(primarySofter);
      double primaryDeltaR = primaryHarder.delta_R(primarySofter);
      double primaryZ =
          primarySofter.modp() / (primarySofter.modp() + primaryHarder.modp());
      double primaryKt = primarySofter.modp() * primaryDeltaR;
      // PseudoJet primaryHarder = (j1.pt() > j2.pt()) ? j1 : j2;
      // PseudoJet primarySofter = (j1.pt() > j2.pt()) ? j2 : j1;
      // double primaryTheta = primaryHarder.delta_phi_to(primarySofter);
      // double primaryDeltaR = primaryHarder.delta_R(primarySofter);
      // double primaryZ =
      //     primarySofter.modp() / (primarySofter.modp() + primaryHarder.modp());
      // double primaryKt = primarySofter.pt() * primaryDeltaR;
      // double primaryKt = primarySofter.modp() * sin(primaryTheta);
      // Update maxPrimary if this kt is larger
      bool isselect = false;
      std::vector<PseudoJet> vec1, vec2, vec3, vec4;
      if (!primaryHarder_harder.has_structure())
        vec1 = {};
      else
        vec1 = primaryHarder_harder.constituents();
      if (!primaryHarder_softer.has_structure())
        vec2 = {};
      else
        vec2 = primaryHarder_softer.constituents();
      if (!primarySofter_harder.has_structure())
        vec3 = {};
      else
        vec3 = primarySofter_harder.constituents();
      if (!primarySofter_softer.has_structure())
        vec4 = {};
      else
        vec4 = primarySofter_softer.constituents();
      if (primaryKt > maxPrimary.kt && primaryZ > z1cut && primaryKt < kt1cut)
      {
        isselect = true;
        maxPrimary = {primaryHarder,
                      primarySofter,
                      primaryHarder_harder,
                      primaryHarder_softer,
                      primarySofter_harder,
                      primarySofter_softer,
                      PseudoJetToTLorentzVector(primaryHarder),
                      PseudoJetToTLorentzVector(primarySofter),
                      primaryTheta,
                      primaryDeltaR,
                      primaryZ,
                      primaryKt,
                      TVector3(),
                      primaryHarder.constituents(),
                      primarySofter.constituents(),
                      vec1,
                      vec2,
                      vec3,
                      vec4};
        maxSecondary = {
            PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), TLorentzVector(), TLorentzVector(), 0.0, 0.0,
            0.0, 0.0, TVector3()};
        maxThirdary = {
            PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), TLorentzVector(), TLorentzVector(), 0.0, 0.0,
            0.0, 0.0, TVector3()};
        PseudoJet secondplanejet = primaryHarder;
        PseudoJet thirdplanejet = primaryHarder;
        if (issecondsoft)
          secondplanejet = primarySofter;
        while (secondplanejet.has_parents(j3, j4))
        {
          PseudoJet secondaryHarder = (j3.modp() > j4.modp()) ? j3 : j4;
          PseudoJet secondarySofter = (j3.modp() > j4.modp()) ? j4 : j3;
          PseudoJet secondaryHarder_harder, secondaryHarder_softer, j31, j32;
          if (secondaryHarder.has_parents(j31, j32))
          {
            secondaryHarder_harder = (j31.modp() > j32.modp()) ? j31 : j32;
            secondaryHarder_softer = (j31.modp() > j32.modp()) ? j32 : j31;
          }
          PseudoJet secondarySofter_harder, secondarySofter_softer, j41, j42;
          if (secondarySofter.has_parents(j41, j42))
          {
            secondarySofter_harder = (j41.modp() > j42.modp()) ? j41 : j42;
            secondarySofter_softer = (j41.modp() > j42.modp()) ? j42 : j41;
          }
          double secondaryTheta = secondaryHarder.delta_phi_to(secondarySofter);
          double secondaryDeltaR =
              secondaryHarder.delta_R(secondarySofter);
          double secondaryZ = secondarySofter.modp() /
                              (secondarySofter.modp() + secondaryHarder.modp());
          double secondaryKt = secondarySofter.modp() * secondaryDeltaR;
          std::vector<PseudoJet> vec1, vec2, vec3, vec4;
          if (!secondaryHarder_harder.has_structure())
            vec1 = {};
          else
            vec1 = secondaryHarder_harder.constituents();
          if (!secondaryHarder_softer.has_structure())
            vec2 = {};
          else
            vec2 = secondaryHarder_softer.constituents();
          if (!secondarySofter_harder.has_structure())
            vec3 = {};
          else
            vec3 = secondarySofter_harder.constituents();
          if (!secondarySofter_softer.has_structure())
            vec4 = {};
          else
            vec4 = secondarySofter_softer.constituents();
          if (secondaryKt > maxSecondary.kt && secondaryZ > z2cut && secondaryKt > kt2cut)
          {
            maxSecondary = {secondaryHarder,
                            secondarySofter,
                            secondaryHarder_harder,
                            secondaryHarder_softer,
                            secondarySofter_harder,
                            secondarySofter_softer,
                            PseudoJetToTLorentzVector(secondaryHarder),
                            PseudoJetToTLorentzVector(secondarySofter),
                            secondaryTheta,
                            secondaryDeltaR,
                            secondaryZ,
                            secondaryKt,
                            TVector3(),
                            secondaryHarder.constituents(),
                            secondarySofter.constituents(),
                            vec1,
                            vec2,
                            vec3,
                            vec4};
          }
          secondplanejet = secondaryHarder;
        }
        while (thirdplanejet.has_parents(j3, j4))
        {
          PseudoJet thirdaryHarder = (j3.modp() > j4.modp()) ? j3 : j4;
          PseudoJet thirdarySofter = (j3.modp() > j4.modp()) ? j4 : j3;
          PseudoJet thirdaryHarder_harder, thirdaryHarder_softer, j31, j32;
          if (thirdaryHarder.has_parents(j31, j32))
          {
            thirdaryHarder_harder = (j31.modp() > j32.modp()) ? j31 : j32;
            thirdaryHarder_softer = (j31.modp() > j32.modp()) ? j32 : j31;
          }
          PseudoJet thirdarySofter_harder, thirdarySofter_softer, j41, j42;
          if (thirdarySofter.has_parents(j41, j42))
          {
            thirdarySofter_harder = (j41.modp() > j42.modp()) ? j41 : j42;
            thirdarySofter_softer = (j41.modp() > j42.modp()) ? j42 : j41;
          }
          double thirdaryTheta = thirdaryHarder.delta_phi_to(thirdarySofter);
          double thirdaryDeltaR =
              thirdaryHarder.delta_R(thirdarySofter);
          double thirdaryZ = thirdarySofter.modp() /
                             (thirdarySofter.modp() + thirdaryHarder.modp());
          double thirdaryKt = thirdarySofter.modp() * thirdaryDeltaR;

          // PseudoJet thirdaryHarder = (j3.pt() > j4.pt()) ? j3 : j4;
          // PseudoJet thirdarySofter = (j3.pt() > j4.pt()) ? j4 : j3;
          // double thirdaryTheta = thirdaryHarder.delta_phi_to(thirdarySofter);
          // double thirdaryDeltaR =
          //     thirdaryHarder.delta_R(thirdarySofter);
          // double thirdaryZ = thirdarySofter.modp() /
          //                    (thirdarySofter.modp() + thirdaryHarder.modp());
          // double thirdaryKt = thirdarySofter.pt() * thirdaryDeltaR;
          std::vector<PseudoJet> vec1, vec2, vec3, vec4;
          if (!thirdaryHarder_harder.has_structure())
            vec1 = {};
          else
            vec1 = thirdaryHarder_harder.constituents();
          if (!thirdaryHarder_softer.has_structure())
            vec2 = {};
          else
            vec2 = thirdaryHarder_softer.constituents();
          if (!thirdarySofter_harder.has_structure())
            vec3 = {};
          else
            vec3 = thirdarySofter_harder.constituents();
          if (!thirdarySofter_softer.has_structure())
            vec4 = {};
          else
            vec4 = thirdarySofter_softer.constituents();
          if (thirdaryKt > maxThirdary.kt && thirdaryZ > z2cut && thirdaryKt > kt2cut)
          {
            maxThirdary = {thirdaryHarder,
                           thirdarySofter,
                           thirdaryHarder_harder,
                           thirdaryHarder_softer,
                           thirdarySofter_harder,
                           thirdarySofter_softer,
                           PseudoJetToTLorentzVector(thirdaryHarder),
                           PseudoJetToTLorentzVector(thirdarySofter),
                           thirdaryTheta,
                           thirdaryDeltaR,
                           thirdaryZ,
                           thirdaryKt,
                           TVector3(),
                           thirdaryHarder.constituents(),
                           thirdarySofter.constituents(),
                           vec1,
                           vec2,
                           vec3,
                           vec4};
          }
          thirdplanejet = thirdaryHarder;
        }
      }
      isselects.push_back(isselect);
      j0 = primaryHarder; // Move to the next iteration with the harder subjet
    }

    int index = -1;
    for (int index0 = isselects.size() - 1; index0 >= 0; index0--)
    {
      if (isselects.at(index0) == true)
      {
        index = index0 + 1;
        break;
      }
    }
    j0 = jinit;
    JetBranch::BChadronPlanes bcplanes = findBhadronPlanes(j0, particlesinfo);
    /// Add plane vector
    maxPrimary.GetVector(particlesinfo);
    maxSecondary.GetVector(particlesinfo);
    maxThirdary.GetVector(particlesinfo);
    threeplanes threeplanes_out;
    maxPrimary.initJet.SetPtEtaPhiE(jinit.pt(), jinit.eta(), jinit.phi(), jinit.e());
    maxSecondary.initJet.SetPtEtaPhiE(jinit.pt(), jinit.eta(), jinit.phi(), jinit.e());
    maxThirdary.initJet.SetPtEtaPhiE(jinit.pt(), jinit.eta(), jinit.phi(), jinit.e());
    threeplanes_out.first = maxPrimary;
    threeplanes_out.second = maxSecondary;
    threeplanes_out.third = maxThirdary;
    threeplanes_out.deltaPhi = JetBranch::DeltaPhi(maxPrimary, maxSecondary);
    threeplanes_out.deltaPhi2 = JetBranch::DeltaPhi(maxPrimary, maxThirdary);
    threeplanes_out.firstindex = index;
    threeplanes_out.bchadrons = bcplanes;
    return threeplanes_out;
  }
  static int netflavour(PseudoJet j, vector<ParticleInfo> particlesinfo, double minflavourpt = 0)
  {
    vector<int> pdgids;
    if (j == 0)
      return -1;
    auto particles = j.constituents();
    for (int i = 0; i < particles.size(); i++)
    {
      PseudoJet particle = particles.at(i);
      if (particle.pt() < minflavourpt)
        continue;
      pdgids.push_back(particlesinfo.at(particle.user_index()).pdgid);
    }
    pdgids.erase(std::remove(pdgids.begin(), pdgids.end(), 21), pdgids.end());
    std::unordered_map<int, int> count_map;
    for (int num : pdgids)
    {
      count_map[num]++;
    }

    pdgids.clear();
    for (auto &element : count_map)
    {
      int num = element.first;
      int count = element.second;
      if (count > 0)
      {
        int opposite = -num;
        if (count_map.count(opposite) && count_map[opposite] > 0)
        {
          int minCount = std::min(count, count_map[opposite]);
          count_map[num] -= minCount;
          count_map[opposite] -= minCount;
        }

        // 将剩余的未抵消的元素添加回vec中
        while (count_map[num] > 0)
        {
          pdgids.push_back(num);
          count_map[num]--;
        }
      }
    }

    if (pdgids.size() == 0)
      return 21;
    if (pdgids.size() == 1)
      return pdgids.at(0);
    if (pdgids.size() > 1)
      return 0;
    return 0;
  }
  static bool isqq(int pdgid1, int pdgid2, int pdgid3 = -1)
  {
    if (pdgid1 != 0 && pdgid1 + pdgid2 == 0 && pdgid3 == -1)
      return true;
    if (pdgid1 != 0 && pdgid1 + pdgid2 == 0 && pdgid3 == 21)
      return true;
    return false;
  }
  static bool isgg(int pdgid1, int pdgid2, int pdgid3 = -1)
  {
    if (pdgid1 == 21 && pdgid2 == 21 && pdgid3 == -1)
      return true;
    if (pdgid1 == 21 && pdgid2 == 21 && pdgid3 == 21)
      return true;
    return false;
  }
  static bool matchPlanes(JetBranch::PlaneVariables &src_plane,
                          JetBranch::PlaneVariables &dst_plane, double drcutweight = 1)
  {

    std::vector<TLorentzVector> src = {src_plane.harder_vect,
                                       src_plane.softer_vect};
    std::vector<TLorentzVector> dst = {dst_plane.harder_vect,
                                       dst_plane.softer_vect};

    double drcut = src[0].DeltaR(src[1]) * drcutweight;
    // match_result[src] = dst
    auto match_result = Matching::matchJet(src, dst, drcut);
    bool allmatch = std::all_of(match_result.begin(), match_result.end(),
                                [](int i)
                                { return i > -1; });

    if (allmatch)
    {
      src_plane.isqq = dst_plane.isqq;
      src_plane.isgg = dst_plane.isgg;
      src_plane.softer_flav = dst_plane.softer_flav;
      src_plane.harder_flav = dst_plane.harder_flav;

      src_plane.harder_init_descri = dst_plane.harder_init_descri;
      src_plane.harder_init_pdgid = dst_plane.harder_init_pdgid;
      src_plane.harder_final_descri = dst_plane.harder_final_descri;
      src_plane.harder_final_pdgid = dst_plane.harder_final_pdgid;
      src_plane.softer_init_descri = dst_plane.softer_init_descri;
      src_plane.softer_init_pdgid = dst_plane.softer_init_pdgid;
      src_plane.softer_final_descri = dst_plane.softer_final_descri;
      src_plane.softer_final_pdgid = dst_plane.softer_final_pdgid;
    }
    return allmatch;
  }

  static vector<std::vector<int>>
  matchPlanes(vector<std::vector<JetBranch::twoplanes>> &src_twoplanes,
              vector<std::vector<JetBranch::twoplanes>> &dst_twoplanes,
              TString opt = "second", double drcutweight = 1)
  {
    std::vector<JetBranch::PlaneVariables> src_plane, dst_plane;
    vector<std::vector<int>> match_results;
    if (opt.EqualTo("second"))
    {
      for (size_t i = 0; i < src_twoplanes.size(); i++)
      {
        std::vector<int> match_results0;
        for (size_t j = 0; j < src_twoplanes.at(i).size(); j++)
        {
          match_results0.push_back(-1);
          src_twoplanes[i][j].second.isgg = false;
          src_twoplanes[i][j].second.isqq = false;
          src_twoplanes[i][j].first.isgg = false;
          src_twoplanes[i][j].first.isqq = false;
          for (size_t k = 0; k < dst_twoplanes.at(i).size(); k++)
          {
            if (matchPlanes(src_twoplanes[i][j].second, dst_twoplanes[i][k].second, drcutweight))
            {
              match_results0[j] = k;
            }
          }
        }
        match_results.push_back(match_results0);
      }
    }
    if (opt.EqualTo("full"))
    {
      for (size_t i = 0; i < src_twoplanes.size(); i++)
      {
        std::vector<int> match_results0;
        for (size_t j = 0; j < src_twoplanes.at(i).size(); j++)
        {
          match_results0.push_back(-1);
          for (size_t k = 0; k < dst_twoplanes.at(i).size(); k++)
          {
            if (matchPlanes(src_twoplanes[i][j].first, dst_twoplanes[i][k].first, drcutweight) && matchPlanes(src_twoplanes[i][j].second, dst_twoplanes[i][k].second, drcutweight))
            {
              match_results0[j] = k;
            }
          }
        }
        match_results.push_back(match_results0);
      }
    }
    return match_results;
  }
  static std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
  matchPlanes(vector<std::vector<JetBranch::threeplanes>> &src_threeplanes,
              vector<std::vector<JetBranch::threeplanes>> &dst_threeplanes,
              TString opt = "second", double drcutweight = 1)
  {
    std::vector<JetBranch::PlaneVariables> src_plane, dst_plane;
    vector<std::vector<int>> match_results2;
    vector<std::vector<int>> match_results3;
    if (opt.EqualTo("second"))
    {
      for (size_t i = 0; i < src_threeplanes.size(); i++)
      {
        std::vector<int> match_results02;
        std::vector<int> match_results03;
        for (size_t j = 0; j < src_threeplanes.at(i).size(); j++)
        {
          match_results02.push_back(-1);
          match_results03.push_back(-1);
          src_threeplanes[i][j].first.isgg = false;
          src_threeplanes[i][j].first.isqq = false;
          src_threeplanes[i][j].second.isgg = false;
          src_threeplanes[i][j].second.isqq = false;
          src_threeplanes[i][j].third.isgg = false;
          src_threeplanes[i][j].third.isqq = false;

          src_threeplanes[i][j].first.harder_flav = 0;
          src_threeplanes[i][j].first.softer_flav = 0;
          src_threeplanes[i][j].second.harder_flav = 0;
          src_threeplanes[i][j].second.softer_flav = 0;
          src_threeplanes[i][j].third.harder_flav = 0;
          src_threeplanes[i][j].third.softer_flav = 0;

          src_threeplanes[i][j].first.harder_init_pdgid = 0;
          src_threeplanes[i][j].first.harder_init_descri = 0;
          src_threeplanes[i][j].first.harder_final_pdgid = 0;
          src_threeplanes[i][j].first.harder_final_descri = 0;
          src_threeplanes[i][j].first.softer_init_pdgid = 0;
          src_threeplanes[i][j].first.softer_init_descri = 0;
          src_threeplanes[i][j].first.softer_final_pdgid = 0;
          src_threeplanes[i][j].first.softer_final_descri = 0;

          src_threeplanes[i][j].second.harder_init_pdgid = 0;
          src_threeplanes[i][j].second.harder_init_descri = 0;
          src_threeplanes[i][j].second.harder_final_pdgid = 0;
          src_threeplanes[i][j].second.harder_final_descri = 0;
          src_threeplanes[i][j].second.softer_init_pdgid = 0;
          src_threeplanes[i][j].second.softer_init_descri = 0;
          src_threeplanes[i][j].second.softer_final_pdgid = 0;
          src_threeplanes[i][j].second.softer_final_descri = 0;

          src_threeplanes[i][j].third.harder_init_pdgid = 0;
          src_threeplanes[i][j].third.harder_init_descri = 0;
          src_threeplanes[i][j].third.harder_final_pdgid = 0;
          src_threeplanes[i][j].third.harder_final_descri = 0;
          src_threeplanes[i][j].third.softer_init_pdgid = 0;
          src_threeplanes[i][j].third.softer_init_descri = 0;
          src_threeplanes[i][j].third.softer_final_pdgid = 0;
          src_threeplanes[i][j].third.softer_final_descri = 0;

          for (size_t k = 0; k < dst_threeplanes.at(i).size(); k++)
          {
            matchPlanes(src_threeplanes[i][j].first, dst_threeplanes[i][k].first, drcutweight);
            if (matchPlanes(src_threeplanes[i][j].second, dst_threeplanes[i][k].second, drcutweight) ||
                matchPlanes(src_threeplanes[i][j].second, dst_threeplanes[i][k].third, drcutweight))
            {
              match_results02[j] = k;
            }
            if (matchPlanes(src_threeplanes[i][j].third, dst_threeplanes[i][k].third, drcutweight) ||
                matchPlanes(src_threeplanes[i][j].third, dst_threeplanes[i][k].second, drcutweight))
            {
              match_results03[j] = k;
            }
          }
        }
        match_results2.push_back(match_results02);
        match_results3.push_back(match_results03);
      }
    }
    return std::make_pair(match_results2, match_results3);
  }
  static vector<std::vector<JetBranch::threeplanes>>
  GetmatchloopPlanes(vector<std::vector<JetBranch::threeplanes>> &src_threeplanes,
                     vector<std::vector<std::vector<ParticleInfo>>> &dst_pseudojets,
                     TString opt = "second", double drcutweight = 1)
  {
    std::vector<std::vector<JetBranch::threeplanes>> dst_threeplanes = src_threeplanes;
    for (auto &sub_vector : dst_threeplanes)
    {
      sub_vector.clear(); // 清空每个子 vector
    }
    std::vector<JetBranch::PlaneVariables> src_plane, dst_plane;
    vector<std::vector<int>> match_results2;
    vector<std::vector<int>> match_results3;
    if (opt.EqualTo("second"))
    {
      for (size_t i = 0; i < src_threeplanes.size(); i++)
      {
        std::vector<int> match_results02;
        std::vector<int> match_results03;
        for (size_t j = 0; j < src_threeplanes.at(i).size(); j++)
        {
          match_results02.push_back(-1);
          match_results03.push_back(-1);
          src_threeplanes[i][j].first.isgg = false;
          src_threeplanes[i][j].first.isqq = false;
          src_threeplanes[i][j].second.isgg = false;
          src_threeplanes[i][j].second.isqq = false;
          src_threeplanes[i][j].third.isgg = false;
          src_threeplanes[i][j].third.isqq = false;

          src_threeplanes[i][j].first.harder_flav = 0;
          src_threeplanes[i][j].first.softer_flav = 0;
          src_threeplanes[i][j].second.harder_flav = 0;
          src_threeplanes[i][j].second.softer_flav = 0;
          src_threeplanes[i][j].third.harder_flav = 0;
          src_threeplanes[i][j].third.softer_flav = 0;

          src_threeplanes[i][j].first.harder_init_pdgid = 0;
          src_threeplanes[i][j].first.harder_init_descri = 0;
          src_threeplanes[i][j].first.harder_final_pdgid = 0;
          src_threeplanes[i][j].first.harder_final_descri = 0;
          src_threeplanes[i][j].first.softer_init_pdgid = 0;
          src_threeplanes[i][j].first.softer_init_descri = 0;
          src_threeplanes[i][j].first.softer_final_pdgid = 0;
          src_threeplanes[i][j].first.softer_final_descri = 0;

          src_threeplanes[i][j].second.harder_init_pdgid = 0;
          src_threeplanes[i][j].second.harder_init_descri = 0;
          src_threeplanes[i][j].second.harder_final_pdgid = 0;
          src_threeplanes[i][j].second.harder_final_descri = 0;
          src_threeplanes[i][j].second.softer_init_pdgid = 0;
          src_threeplanes[i][j].second.softer_init_descri = 0;
          src_threeplanes[i][j].second.softer_final_pdgid = 0;
          src_threeplanes[i][j].second.softer_final_descri = 0;

          src_threeplanes[i][j].third.harder_init_pdgid = 0;
          src_threeplanes[i][j].third.harder_init_descri = 0;
          src_threeplanes[i][j].third.harder_final_pdgid = 0;
          src_threeplanes[i][j].third.harder_final_descri = 0;
          src_threeplanes[i][j].third.softer_init_pdgid = 0;
          src_threeplanes[i][j].third.softer_init_descri = 0;
          src_threeplanes[i][j].third.softer_final_pdgid = 0;
          src_threeplanes[i][j].third.softer_final_descri = 0;

          for (size_t k = 0; k < dst_pseudojets.at(i).size(); k++)
          {
            std::vector<ParticleInfo> constituents = dst_pseudojets.at(i).at(k);
            double z1cut = 0;
            double z2cut = 0;
            double kt1cut = 1;
            double kt2cut = 0;
            double issecondsoft = true;
            vector<JetBranch::twoplanes> twoplanes;
            std::vector<PseudoJet> particles;
            std::vector<ParticleInfo> particlesinfo;
            for (size_t cs = 0; cs < constituents.size(); ++cs)
            {
              TLorentzVector p;
              p.SetPtEtaPhiE(constituents.at(cs).pt, constituents.at(cs).eta,
                             constituents.at(cs).phi, constituents.at(cs).e);
              PseudoJet particle = PseudoJet(p.Px(), p.Py(), p.Pz(), p.Energy());
              int pdgid = constituents.at(i).pdgid;
              ParticleInfo particleInfo(pdgid, constituents.at(i).chargeInt);
              particlesinfo.push_back(particleInfo);
              particle.set_user_index(cs);
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
            PseudoJet j0 = jets_IFN[0];
            PseudoJet j1, j2;
            std::vector<std::pair<PseudoJet, PseudoJet>> looppseudojets;
            if (j0.has_parents(j1, j2))
              looppseudojets = {{j1, j2}};
            while (looppseudojets.size() != 0)
            {
              for (int pseudo = 0; pseudo < looppseudojets.size(); pseudo++)
              {
                PseudoJet harder = (looppseudojets.at(pseudo).first.modp() > looppseudojets.at(pseudo).second.modp()) ? looppseudojets.at(pseudo).first : looppseudojets.at(pseudo).second;
                PseudoJet softer = (looppseudojets.at(pseudo).first.modp() > looppseudojets.at(pseudo).second.modp()) ? looppseudojets.at(pseudo).second : looppseudojets.at(pseudo).first;
                double deltaR = PseudoJetToTLorentzVector(harder).DeltaR(PseudoJetToTLorentzVector(softer));
                int pdgid1 = JetBranch::GetIFNFlavour(harder);
                int pdgid2 = JetBranch::GetIFNFlavour(softer);
                // clang-format off
                PlaneVariables pseudoplanes = {
                    harder, softer, PseudoJet(), PseudoJet(), PseudoJet(), PseudoJet(), 
                    PseudoJetToTLorentzVector(harder), PseudoJetToTLorentzVector(softer), 
                    0, deltaR, 0, 0, {0, 0, 0}, harder.constituents(), softer.constituents(), {}, {}, {}, {},
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {}, {}, {}, {}, {}, {}, 
                    GetIFNFlavour(harder), GetIFNFlavour(softer), 
                    GetIFNPdgid(softer, false), GetIFNDescri(softer, false), 
                    GetIFNPdgid(softer, true), GetIFNDescri(softer, true), 
                    GetIFNPdgid(harder, false), GetIFNDescri(harder, false), 
                    GetIFNPdgid(harder, true), GetIFNDescri(harder, true), 
                    0, 0, 0, 0};
                // clang-format on
                if (matchPlanes(src_threeplanes[i][j].second, pseudoplanes, drcutweight))
                {
                  JetBranch::threeplanes pseudo_threeplanes;
                  pseudo_threeplanes.second = pseudoplanes;
                  dst_threeplanes.at(i).push_back(pseudo_threeplanes);
                  break;
                }
              }
              if (dst_threeplanes.at(i).size() != 0)
                break;
              std::vector<std::pair<PseudoJet, PseudoJet>> temp;
              for (int pseudo = 0; pseudo < looppseudojets.size(); pseudo++)
              {
                PseudoJet j3, j4, j5, j6;
                if (looppseudojets.at(pseudo).first.has_parents(j3, j4))
                {
                  temp.push_back({j3, j4});
                }
                if (looppseudojets.at(pseudo).second.has_parents(j5, j6))
                {
                  temp.push_back({j5, j6});
                }
              }
              looppseudojets = temp;
            }
          }
        }
      }
    }
    return dst_threeplanes;
  }
  static std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
  matchPlanesLeadingJet(vector<std::vector<JetBranch::threeplanes>> &src_threeplanes,
                        vector<std::vector<JetBranch::threeplanes>> &dst_threeplanes,
                        TString opt = "second", double drcutweight = 1)
  {
    std::vector<JetBranch::PlaneVariables> src_plane, dst_plane;
    vector<std::vector<int>> match_results2;
    vector<std::vector<int>> match_results3;
    if (opt.EqualTo("second"))
    {
      for (size_t i = 0; i < src_threeplanes.size(); i++)
      {
        std::vector<int> match_results02;
        std::vector<int> match_results03;
        for (size_t j = 0; j < src_threeplanes.at(i).size(); j++)
        {
          match_results02.push_back(-1);
          match_results03.push_back(-1);
          src_threeplanes[i][j].first.isgg = false;
          src_threeplanes[i][j].first.isqq = false;
          src_threeplanes[i][j].second.isgg = false;
          src_threeplanes[i][j].second.isqq = false;
          src_threeplanes[i][j].third.isgg = false;
          src_threeplanes[i][j].third.isqq = false;

          src_threeplanes[i][j].first.harder_flav = 0;
          src_threeplanes[i][j].first.softer_flav = 0;
          src_threeplanes[i][j].second.harder_flav = 0;
          src_threeplanes[i][j].second.softer_flav = 0;
          src_threeplanes[i][j].third.harder_flav = 0;
          src_threeplanes[i][j].third.softer_flav = 0;

          src_threeplanes[i][j].first.harder_init_pdgid = 0;
          src_threeplanes[i][j].first.harder_init_descri = 0;
          src_threeplanes[i][j].first.harder_final_pdgid = 0;
          src_threeplanes[i][j].first.harder_final_descri = 0;
          src_threeplanes[i][j].first.softer_init_pdgid = 0;
          src_threeplanes[i][j].first.softer_init_descri = 0;
          src_threeplanes[i][j].first.softer_final_pdgid = 0;
          src_threeplanes[i][j].first.softer_final_descri = 0;

          src_threeplanes[i][j].second.harder_init_pdgid = 0;
          src_threeplanes[i][j].second.harder_init_descri = 0;
          src_threeplanes[i][j].second.harder_final_pdgid = 0;
          src_threeplanes[i][j].second.harder_final_descri = 0;
          src_threeplanes[i][j].second.softer_init_pdgid = 0;
          src_threeplanes[i][j].second.softer_init_descri = 0;
          src_threeplanes[i][j].second.softer_final_pdgid = 0;
          src_threeplanes[i][j].second.softer_final_descri = 0;

          src_threeplanes[i][j].third.harder_init_pdgid = 0;
          src_threeplanes[i][j].third.harder_init_descri = 0;
          src_threeplanes[i][j].third.harder_final_pdgid = 0;
          src_threeplanes[i][j].third.harder_final_descri = 0;
          src_threeplanes[i][j].third.softer_init_pdgid = 0;
          src_threeplanes[i][j].third.softer_init_descri = 0;
          src_threeplanes[i][j].third.softer_final_pdgid = 0;
          src_threeplanes[i][j].third.softer_final_descri = 0;

          for (size_t k = 0; k < dst_threeplanes.at(i).size(); k++)
          {
            matchPlanes(src_threeplanes[i][j].first, dst_threeplanes[i][k].first, drcutweight);
            if (matchPlanes(src_threeplanes[i][j].second, dst_threeplanes[i][k].second, drcutweight) ||
                matchPlanes(src_threeplanes[i][j].second, dst_threeplanes[i][k].third, drcutweight))
            {
              match_results02[j] = k;
            }
            if (matchPlanes(src_threeplanes[i][j].third, dst_threeplanes[i][k].third, drcutweight) ||
                matchPlanes(src_threeplanes[i][j].third, dst_threeplanes[i][k].second, drcutweight))
            {
              match_results03[j] = k;
            }
          }
        }
        match_results2.push_back(match_results02);
        match_results3.push_back(match_results03);
      }
    }
    return std::make_pair(match_results2, match_results3);
  }
  static void SetUserIndex(PseudoJet &jet, int &index)
  {
    if (jet.user_index() != -1)
    {
      // 如果已经设置了index，就不再处理
      return;
    }
    jet.set_user_index(index++); // 设置当前jet的index，并增加index

    PseudoJet parent1, parent2;
    if (jet.has_parents(parent1, parent2))
    {
      SetUserIndex(parent1, index); // 递归标记第一个parent
      SetUserIndex(parent2, index); // 递归标记第二个parent
    }
  }

  static void SetPseudoCAJetsIndex(vector<PseudoJet> &jets)
  {
    int k = 0;
    for (PseudoJet &jet : jets)
    {
      SetUserIndex(jet, k);
    }
  }
  static int GetIFNFlavour(PseudoJet &jet)
  {
    if (jet.pt() == 0)
      return 0;
    int pdgid = FlavHistory::current_flavour_of(jet).pdg_code();
    // std::cout << pdgid << std::endl;
    if (pdgid == 0)
    {
      TString flav = FlavHistory::current_flavour_of(jet).description();
      if (flav.EqualTo("[g]"))
      {
        pdgid = 21;
      }
    }
    return pdgid;
  }

  static int GetNtracks(vector<PseudoJet> constituents, vector<ParticleInfo> particlesinfo)
  {
    int ncharge = 0;
    int x = 0;
    int y = 0;
    for (auto &particle : constituents)
    {
      int particle_index = particle.user_index();

      double charge = particlesinfo.at(particle_index).charge;
      if (abs(charge) >= 0.1)
        ncharge++;
    }
    return ncharge;
  }

  static int GetIFNPdgid(PseudoJet &jet, bool iscurrent = true)
  {
    if (jet.pt() == 0)
      return 0;
    int pdgid = 0;
    if (iscurrent)
      pdgid = FlavHistory::current_flavour_of(jet).pdg_code();
    else
      pdgid = FlavHistory::initial_flavour_of(jet).pdg_code();
    // std::cout << pdgid << std::endl;
    return pdgid;
  };

  static int GetIFNDescri(PseudoJet &jet, bool iscurrent = true)
  {
    if (jet.pt() == 0)
      return 0;
    int pdgid;
    TString flav;
    if (iscurrent)
      flav = FlavHistory::current_flavour_of(jet).description();
    else
      flav = FlavHistory::initial_flavour_of(jet).description();

    if (flav.EqualTo("[g]"))
      pdgid = 21;
    else if (flav.EqualTo("[d ]"))
      pdgid = 1;
    else if (flav.EqualTo("[u ]"))
      pdgid = 2;
    else if (flav.EqualTo("[s ]"))
      pdgid = 3;
    else if (flav.EqualTo("[c ]"))
      pdgid = 4;
    else if (flav.EqualTo("[b ]"))
      pdgid = 5;
    else if (flav.EqualTo("[t ]"))
      pdgid = 6;
    else if (flav.EqualTo("[dbar ]"))
      pdgid = -1;
    else if (flav.EqualTo("[ubar ]"))
      pdgid = -2;
    else if (flav.EqualTo("[sbar ]"))
      pdgid = -3;
    else if (flav.EqualTo("[cbar ]"))
      pdgid = -4;
    else if (flav.EqualTo("[bbar ]"))
      pdgid = -5;
    else if (flav.EqualTo("[tbar ]"))
      pdgid = -6;
    else
      pdgid = 0;

    return pdgid;
  };

  static int GetGhostNumber(PseudoJet &j0, std::vector<ParticleInfo> particlesinfo, int type = 5)
  {
    int bnum = 0;
    int cnum = 0;
    PseudoJet jinit = j0;
    for (int i = 0; i < jinit.constituents().size(); i++)
    {
      auto particle = jinit.constituents().at(i);
      if ((abs(particlesinfo.at(particle.user_index()).pdgid / 100 % 10) == 5 ||
           abs(particlesinfo.at(particle.user_index()).pdgid / 1000 % 10) == 5))
        bnum++;
      else if ((abs(particlesinfo.at(particle.user_index()).pdgid / 100 % 10) == 4 ||
                abs(particlesinfo.at(particle.user_index()).pdgid / 1000 % 10) == 4))
        cnum++;
    }
    if (type == 5)
      return bnum;
    else if (type == 4)
      return cnum;
    else
      return 0;
  };
};