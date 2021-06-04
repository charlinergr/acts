// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/GraphCreation/ModuleMapCreatorWriter.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Io/Csv/CsvHitReader.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include <Acts/Geometry/TrackingVolume.hpp>
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <boost/bimap.hpp>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <limits>

using namespace ActsExamples;

void ModuleMapCreatorWriter::examineLayer(const Acts::TrackingVolume* vol) {
  
  // confined layers
  const Acts::LayerArray* clayers = vol->confinedLayers();
  if (clayers) {
    auto objects = clayers->arrayObjects();  // LayerArray
    for (size_t i = 0; i<objects.size(); i++) {
      auto sa = objects[i].get()->surfaceArray();
      if (sa) {
	      auto surfaces = sa->surfaces();
	      for (size_t ii = 0; ii<surfaces.size(); ii++) {
	        if (objects[i].get()->layerType()==Acts::LayerType::active) {	 
            m_ModuleMap.insert({surfaces[ii]->geometryId().value(),m_ModuleMap.size()});
          }
	      }
      }
    }
  }

  // confined volumes
  auto confinedVolumes = vol->confinedVolumes();
  if (confinedVolumes) {
    auto objects = confinedVolumes->arrayObjects();
    for(size_t i = 0; i<objects.size(); i++) {
      examineLayer(objects[i].get());
    }
  }  
  
}// end void examineLayer

ModuleMapCreatorWriter::ModuleMapCreatorWriter(
				     const ModuleMapCreatorWriter::Config& cfg, Acts::Logging::Level lvl)
  : m_cfg(cfg),
    m_world(nullptr),
    m_logger(Acts::getDefaultLogger("ModuleMapCreatorWriter", lvl))

{
  std::cout<<"#############"<<"\n"
         <<"Input options: "<<"\n"
         <<"  #cut on pt = "<<m_cfg.minPt<<"\n"
         <<"  #cut on nhits = "<<m_cfg.minNHits<<"\n"
         <<" with cuts ? "<<m_cfg.giveCutsValues<<"\n"
         <<std::endl;

  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  
  m_world = m_cfg.trackingGeometry->highestTrackingVolume();
  if (not m_world) {
    throw std::invalid_argument("Could not identify the world volume");
  }
  
  //Fill m_ModuleMap with all active geometryIdentifier of the detector
  examineLayer(m_world);
  ACTS_VERBOSE("Found" << m_ModuleMap.size() << "modules" << "]");
 
  m_LinkCounters = new unsigned int*[m_ModuleMap.size()];
  //storage of cut values

  m_module_z0_max = new float*[m_ModuleMap.size()];
  m_module_z0_min = new float*[m_ModuleMap.size()];
  m_module_phislope_max = new float*[m_ModuleMap.size()];
  m_module_phislope_min = new float*[m_ModuleMap.size()];
  m_module_deta_max = new float*[m_ModuleMap.size()];
  m_module_deta_min = new float*[m_ModuleMap.size()];
  m_module_dphi_max = new float*[m_ModuleMap.size()];
  m_module_dphi_min = new float*[m_ModuleMap.size()];

  for(size_t i = 0; i < m_ModuleMap.size(); i++) {
    m_LinkCounters[i] = new unsigned int[m_ModuleMap.size()];
    m_module_z0_max[i] = new float [m_ModuleMap.size()];
    m_module_z0_min[i] = new float[m_ModuleMap.size()];
    m_module_phislope_max[i] = new float[m_ModuleMap.size()];
    m_module_phislope_min[i] = new float[m_ModuleMap.size()];
    m_module_deta_max[i] = new float[m_ModuleMap.size()];
    m_module_deta_min[i] = new float[m_ModuleMap.size()];
    m_module_dphi_max[i] = new float[m_ModuleMap.size()];
    m_module_dphi_min[i] = new float[m_ModuleMap.size()];
    //initialise all variables with dumb values
    for(size_t y = 0; y < m_ModuleMap.size(); y++) {
      m_module_z0_max[i][y] = -10000000000.;
      m_module_z0_min[i][y] = 10000000000.;
      m_module_phislope_max[i][y] = -10000000000.;
      m_module_phislope_min[i][y] = 10000000000.;
      m_module_deta_max[i][y] = -10000000000.;
      m_module_deta_min[i][y] = 10000000000.;
      m_module_dphi_max[i][y] = -10000000000.;
      m_module_dphi_min[i][y] = 10000000000.;
    }
  }
  
  //Set up the branch
  m_treeModuleMap->Branch("Module1",&m_tree_module1,"Module1/l");
  m_treeModuleMap->Branch("Module2",&m_tree_module2,"Module2/l");
  m_treeModuleMap->Branch("Occurrence",&m_count,"Occurrence/i");
  m_treeModuleMap->Branch("z0_max",&m_z0_max,"z0_max/F");
  m_treeModuleMap->Branch("z0_min",&m_z0_min,"z0_min/F");
  m_treeModuleMap->Branch("phislope_max",&m_phislope_max,"phislope_max/F");
  m_treeModuleMap->Branch("phislope_min",&m_phislope_min,"phislope_min/F");
  m_treeModuleMap->Branch("deta_max",&m_deta_max,"deta_max/F");
  m_treeModuleMap->Branch("deta_min",&m_deta_min,"deta_min/F");
  m_treeModuleMap->Branch("dphi_max",&m_dphi_max,"dphi_max/F");
  m_treeModuleMap->Branch("dphi_min",&m_dphi_min,"dphi_min/F");
}

std::string ModuleMapCreatorWriter::name() const {
  return "ModuleMapCreatorWriter";
}

ProcessCode ModuleMapCreatorWriter::write(const AlgorithmContext& ctx) {
  
  // Get input collections from the event store
  const auto& inputParticles = ctx.eventStore.get<SimParticleContainer>("particles_initial");

  //using ClusterContainer = GeometryIdMultimap<Acts::PlanarModuleCluster>;
  //const auto& inputClusters =ctx.eventStore.get<ClusterContainer>("clusters");
  
  //using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
  //const auto& hitParticlesMap = ctx.eventStore.get<HitParticlesMap>("hit_particles_map");
  // compute particle_id -> {hit_id...} map from the hit_id -> {particle_id...} map on the fly.
  //const auto& particleHitsMap = invertIndexMultimap(hitParticlesMap);

  //using HitsContainer = std::map<Acts::GeometryIdentifier::Value,struct Hitinformation>;
  //const auto& inputHits =ctx.eventStore.get<HitsContainer>("hits");
  using HitsContainer = std::multimap<ActsFatras::Barcode,ActsExamples::CsvHitReader::Hitinformation>;
  const auto& inputHits =ctx.eventStore.get<HitsContainer>("hitsPID");

  // Loop over all particles
  for (const auto& particle : inputParticles) {
    // cut on pt
    if (particle.transverseMomentum()<m_cfg.minPt) {
        continue;
    }

 
   
    auto hits = makeRange(inputHits.equal_range(particle.particleId()));
   
    if (hits.size() <2) continue; //remove particles leaving only one hit here (or no hit)

    //idk how to sort it in a more clean way
    std::vector<ActsExamples::CsvHitReader::Hitinformation> vect_hits;
    for (auto& i : hits) vect_hits.push_back(i.second);
    sort( vect_hits.begin( ), vect_hits.end( ),
      []( ActsExamples::CsvHitReader::Hitinformation& a, ActsExamples::CsvHitReader::Hitinformation & b ){
      return a.index < b.index;
    } );
    
    //Check if two hits of the same particles are on the same module. If True, remove one of them.
    //this is something to check in data sample, should not be possible
    for (size_t i=0; i<vect_hits.size()-1; i++){
      if (vect_hits[i].geometry_id.value() == vect_hits[i+1].geometry_id.value()) vect_hits.erase(vect_hits.begin()+i+1);
    }
    
    //cut on number of hits
    if (vect_hits.size()<m_cfg.minNHits) continue;

    // Fill the matrix of (n_module*n_module) to count numbers of possible connections
    
    for (size_t i=0; i<vect_hits.size()-1; i++){
      auto it1 = m_ModuleMap.left.find(vect_hits[i].geometry_id.value());
      auto it2 = m_ModuleMap.left.find(vect_hits[i+1].geometry_id.value());
      
      if (!(it1 != m_ModuleMap.left.end() && it2 != m_ModuleMap.left.end())) continue; 
      
      size_t index1 = it1->second;
      size_t index2 = it2->second;
      m_LinkCounters[index1][index2]++;
      
      if(!m_cfg.giveCutsValues) continue;
        
      float x1 = vect_hits[i].x;
      float y1 = vect_hits[i].y;
      float z1 = vect_hits[i].z;

      float x2 = vect_hits[i+1].x;
      float y2 = vect_hits[i+1].y;
      float z2 = vect_hits[i+1].z;
      //cut values
      //phi
      float phi1 = atan2(y1,x1);
      float phi2 = atan2(y2,x2);
      float dphi = phi2 - phi1 ;
      if (dphi>TMath::Pi()) dphi-=2*TMath::Pi();
      if (dphi<-TMath::Pi()) dphi+=2*TMath::Pi();
      if (dphi > m_module_dphi_max[index1][index2]) {
        m_module_dphi_max[index1][index2]= dphi; 
      }
      if (dphi < m_module_dphi_min[index1][index2]) {
        m_module_dphi_min[index1][index2] = dphi; 
      }
      //phi slope
      float r1 = sqrt(x1*x1 + y1*y1);
      float r2 = sqrt(x2*x2 + y2*y2);
      float dr = r2 - r1;
      float phi_slope = dphi / dr;
      if (!(isinf(phi_slope)==1)){
        if (phi_slope > m_module_phislope_max[index1][index2]) {
          m_module_phislope_max[index1][index2] = phi_slope; 
        }
        if (phi_slope < m_module_phislope_min[index1][index2]) {
          m_module_phislope_min[index1][index2]= phi_slope; 
        }
      }
      //z0
      float dz = z2 - z1;
      float z0 = z1 - r1 * dz / dr;
      if (!(isinf(z0)==1)){
        if (z0 > m_module_z0_max[index1][index2]) {
          m_module_z0_max[index1][index2] = z0;
        }
        if (z0 < m_module_z0_min[index1][index2]) {
          m_module_z0_min[index1][index2]= z0; 
        }
      }
      //deta
      float r31 = sqrt(r1*r1 + z1*z1);
      float r32 = sqrt(r2*r2 + z2*z2);
      float theta11 = acos(z1/r31);
      float theta22 = acos(z2/r32);
      float eta1 = -log(tan(theta11/2.));
      float eta2 = -log(tan(theta22/2.));
      float deta = eta1 - eta2;
      if (deta > m_module_deta_max[index1][index2]) {
        m_module_deta_max[index1][index2] = deta; 
      }
      if (deta < m_module_deta_min[index1][index2]) {
        m_module_deta_min[index1][index2] = deta; 
      } 
    }
  }
  return ProcessCode::SUCCESS;
}


ProcessCode ModuleMapCreatorWriter::endRun() {

  //Fill output root
  //std::cout.precision(6);
  for (size_t module1 = 0; module1 < m_ModuleMap.size(); module1++) {
    for(size_t module2 = 0; module2 < m_ModuleMap.size(); module2++) {
      if(m_LinkCounters[module1][module2]>0 && module1!=module2 ) {
        m_tree_module1 = m_ModuleMap.right.find(module1)->second;
        m_tree_module2 = m_ModuleMap.right.find(module2)->second;
        m_count = m_LinkCounters[module1][module2];
        m_z0_max = m_module_z0_max[module1][module2];
        m_z0_min = m_module_z0_min[module1][module2];
        
        //std::cout<<std::fixed<<module1<<" "<<module2<<" "<<m_z0_max<<" " << m_z0_min<<std::endl;
        //std::cout<<m_tree_module1 <<" "<<m_tree_module2<<std::endl;
        m_phislope_max = m_module_phislope_max[module1][module2];
        m_phislope_min = m_module_phislope_min[module1][module2];
        m_deta_max = m_module_deta_max[module1][module2];
        m_deta_min = m_module_deta_min[module1][module2];
        m_dphi_max = m_module_dphi_max[module1][module2];
        m_dphi_min = m_module_dphi_min[module1][module2];
        m_treeModuleMap->Fill();
      }
    }
  }

  auto path_RootFile = joinPaths(m_cfg.outputDir, m_cfg.rootName);
  m_RootFile = TFile::Open(path_RootFile.c_str(), "RECREATE");
  m_RootFile->cd();
  m_treeModuleMap->Write();
  m_RootFile->Close();

  return ProcessCode::SUCCESS;
}
