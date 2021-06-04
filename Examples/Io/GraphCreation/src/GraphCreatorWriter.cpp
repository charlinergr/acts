// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "ActsExamples/Io/GraphCreation/GraphCreatorWriter.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/Io/Csv/CsvHitReader.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <boost/graph/graphml.hpp>
#include <boost/graph/adjacency_list.hpp> // we have to test with #include <boost/graph/adjacency_matrix.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <map>
#include "Acts/Utilities/Logger.hpp"

using namespace ActsExamples;

bool GraphCreatorWriter::compare_float(float x, float y){
  float epsilon=0.0001f;
  if(fabs(x - y) < epsilon) return true; //they are same
  return false; //they are not same
}


GraphCreatorWriter::GraphCreatorWriter(
				       const GraphCreatorWriter::Config& cfg, Acts::Logging::Level lvl)
  : m_cfg(cfg),
    m_logger(Acts::getDefaultLogger("GraphCreatorWriter", lvl))

{
std::cout.precision(6);
std::cout<<"#############"<<"\n"
         <<"Input options: "<<"\n"
         <<"  #module map = "<<m_cfg.inputModuleMap<<"\n"
         <<" true graph ? "<<m_cfg.trueGraph<<"\n"
         <<" save graph ? "<<m_cfg.saveGraph<<"\n"
         <<std::endl;

if ( m_cfg.inputModuleMap.empty() ) {
    throw std::invalid_argument("Missing Module Map");
}
}

std::string GraphCreatorWriter::name() const {  
  return "GraphCreatorWriter";
}

ProcessCode GraphCreatorWriter::write(const AlgorithmContext& ctx) {
 
  ACTS_INFO("start");
  std::cout.precision(6);
  using HitsContainer = GeometryIdMultimap<ActsExamples::CsvHitReader::Hitinformation>;
  const auto& inputHits =ctx.eventStore.get<HitsContainer>("hits");

  
  Graph G;
  Graph_true G_true;
  
  //Get the root file containing the module map and all cut values
  TFile *root_ModuleMap = TFile::Open(m_cfg.inputModuleMap.c_str(),"READ");
  TTree *tree_ModuleMap = (TTree*)root_ModuleMap->Get("TreeModuleMap");
  //branch variables
  ULong64_t module1, module2;
  unsigned int m_count;
  float z0_min;
  float phislope_min;
  float deta_min;
  float dphi_min;
  float z0_max;
  float phislope_max;
  float deta_max;
  float dphi_max;
  //get branch
  tree_ModuleMap->SetBranchAddress("Module1",&module1);
  tree_ModuleMap->SetBranchAddress("Module2",&module2);
  tree_ModuleMap->SetBranchAddress("Occurrence",&m_count);
  tree_ModuleMap->SetBranchAddress("z0_max",&z0_max);
  tree_ModuleMap->SetBranchAddress("z0_min",&z0_min);
  tree_ModuleMap->SetBranchAddress("phislope_max",&phislope_max);
  tree_ModuleMap->SetBranchAddress("phislope_min",&phislope_min);
  tree_ModuleMap->SetBranchAddress("deta_max",&deta_max);
  tree_ModuleMap->SetBranchAddress("deta_min",&deta_min);
  tree_ModuleMap->SetBranchAddress("dphi_max",&dphi_max);
  tree_ModuleMap->SetBranchAddress("dphi_min",&dphi_min);

  std::string event_id=std::to_string(ctx.eventNumber);
  std::string filename_id="/event"+event_id+"_ID.txt";
	auto path_id = joinPaths(m_cfg.outputDir, filename_id);              
  std::fstream rausIDfile(path_id,std::ios::out); 

  //check purpose
  int true_edges_number = 0;
  //loop over the module map
  mapB hit_to_node;

  for (int connection=0; connection<tree_ModuleMap->GetEntries(); connection++) {
    tree_ModuleMap->GetEntry(connection);

    //TO DO: possibility to cut on volume and or layer to add in config configurations
    auto hits1 = makeRange(inputHits.equal_range(module1));
    auto hits2 = makeRange(inputHits.equal_range(module2));

    if ( (hits1.size()==0) || (hits2.size()==0) ) continue;
   
    for (auto& i : hits1) {

      float r1 = sqrt(i.second.x*i.second.x + i.second.y*i.second.y);
      float phi1 = atan2(i.second.y,i.second.x);
      float r31 = sqrt(r1*r1 + i.second.z*i.second.z);
      float theta11 = acos(i.second.z/r31);
      float eta1 = -log(tan(theta11/2.));

      for (auto& ii : hits2) {
        //applied cuts and compute geometric information
        //for the case where a true edge matched the cut value, a compare_float function must be used

        float r2 = sqrt(ii.second.x*ii.second.x + ii.second.y*ii.second.y);
        float dr = r2 - r1; //can be = 0, in that case it pass all cuts
        float dz = ii.second.z - i.second.z;
        float z0 = i.second.z - r1 * dz / dr; 

        if ((z0<z0_min || z0>z0_max) && (!isinf(z0)==1)) continue;

        float phi2 = atan2(ii.second.y,ii.second.x);
        float dphi = phi2 - phi1;
        if (dphi>TMath::Pi()) dphi-=2*TMath::Pi();
        if (dphi<-TMath::Pi()) dphi+=2*TMath::Pi();

        if ((dphi<dphi_min || dphi>dphi_max)) continue;

        float phi_slope = dphi / dr;

        if ((phi_slope<phislope_min || phi_slope>phislope_max) &&(!(isinf(phi_slope)==1))) continue;

        float r32 = sqrt(r2*r2 + ii.second.z*ii.second.z);
        float theta22 = acos(ii.second.z/r32);
        float eta2 = -log(tan(theta22/2.));
        float deta = eta1 - eta2;

        if ((deta<deta_min || deta>deta_max)) continue;

        //Turn into graph format
        //Node creation
        
        vertex u=0;
        auto hit1_in_graph = hit_to_node.find(i.second.hit_id);
		    if (hit1_in_graph == hit_to_node.end()){
		      u = boost::add_vertex(G);  
		      G[u].r = r1/1000.;
          G[u].phi =  phi1 /TMath::Pi();
          G[u].z = i.second.z/1000.;
          G[u].hit_id= i.second.hit_id;
          hit_to_node.insert(std::pair<uint64_t,vertex>(i.second.hit_id,u));
          rausIDfile<<i.second.hit_id<<"\n";
        }else{
          u= hit1_in_graph->second;
        }
        
        vertex v=0;
        auto hit2_in_graph = hit_to_node.find(ii.second.hit_id);
		    if (hit2_in_graph == hit_to_node.end()){
		      v = boost::add_vertex(G);  
		      G[v].r = r2/1000.;
          G[v].phi =  phi2/TMath::Pi();
          G[v].z = ii.second.z/1000.;
          G[v].hit_id= ii.second.hit_id;
          hit_to_node.insert(std::pair<uint64_t,vertex>(ii.second.hit_id,v));
          rausIDfile<<ii.second.hit_id<<"\n";
        }else{
          v= hit2_in_graph->second;
        }

        //Edge creation
        edge e; 
        bool b;
        boost::tie(e,b) = boost::add_edge(u,v,G);
        G[e].dEta = deta;
        G[e].dPhi = dphi;    
        G[e].dr = dr;
        G[e].dz= dz;

        if (m_cfg.trueGraph){
          edge_true e_t; 
          boost::tie(e_t,b) = boost::add_edge(u,v,G_true);     
          G_true[u].pt_particle=i.second.pt_particle;
          G_true[u].r=r1;
          G_true[u].PID=i.second.particle_id;
          G_true[u].z=i.second.z;
          G_true[u].eta=eta1;
          G_true[u].phi=phi1;
          G_true[u].hit_id=i.second.hit_id;
          G_true[u].index=i.second.index;

          G_true[v].pt_particle=ii.second.pt_particle;
          G_true[v].r=r2;
          G_true[v].PID=ii.second.particle_id;
          G_true[v].z=ii.second.z;
          G_true[v].eta=eta2;
          G_true[v].phi=phi2;
          G_true[v].hit_id=ii.second.hit_id;
          G_true[v].index=ii.second.index;


          G_true[e_t].is_segment_true=0;
          //Check if the edge is a true one
          if (i.second.particle_id == ii.second.particle_id && i.second.nhits_wihoutDuplicate>=m_cfg.minNHits && i.second.pt_particle>=m_cfg.minPt){
            if (i.second.duplicate == false ){
              if ((i.second.index +1 ==ii.second.index)) {
                G_true[e_t].is_segment_true=1;
                G_true[e_t].pt_particle=ii.second.pt_particle;
                true_edges_number++;
              }
            } else if (i.second.linkUp == true) {
              if ((i.second.index+2 == ii.second.index)){
                G_true[e_t].is_segment_true=1;
                G_true[e_t].pt_particle=ii.second.pt_particle;
                true_edges_number++;
              }
            }
          }
        } //close trueGraph condition
      }
    }
  }
  root_ModuleMap->Close();
  
  std::cout<<ctx.eventNumber<<" "<<"There are "<<boost::num_edges(G)<<" edges"<<std::endl;
  std::cout<<ctx.eventNumber<<" "<<"There are "<<true_edges_number<<" true edges"<<std::endl;
  std::cout<<ctx.eventNumber<<" "<<"There are "<<boost::num_vertices(G)<<" nodes"<<std::endl;
  rausIDfile.close();
  //write in graphML format on request
  
  if (m_cfg.saveGraph){

    std::string evtid=std::to_string(ctx.eventNumber);

    boost::dynamic_properties dp(boost::ignore_other_properties);
    dp.property("r", boost::get(&VertexProperty::r, G));
    dp.property("phi", boost::get(&VertexProperty::phi, G));
    dp.property("z", boost::get(&VertexProperty::z, G));
    dp.property("hit_id", boost::get(&VertexProperty::hit_id, G));
    dp.property("dEta", boost::get(&EdgeProperty::dEta,G));
    dp.property("dPhi", boost::get(&EdgeProperty::dPhi, G));
    dp.property("dr", boost::get(&EdgeProperty::dr, G));
    dp.property("dz", boost::get(&EdgeProperty::dz, G));

    std::string filenameg="/event"+evtid+"_INPUT.txt";
    auto path_graph = joinPaths(m_cfg.outputDir, filenameg);
    std::fstream outGraph(path_graph,std::ios::out); 
    boost::write_graphml( outGraph, G,dp, true);
    outGraph.close();
    if (m_cfg.trueGraph) {
      boost::dynamic_properties dp_t(boost::ignore_other_properties);
      dp_t.property("is_segment_true", boost::get(&EdgePropertyTrue::is_segment_true, G_true));
      dp_t.property("pt_particle", boost::get(&EdgePropertyTrue::pt_particle, G_true));

      dp_t.property("r", boost::get(&VertexPropertyTrue::r, G_true));
      dp_t.property("phi", boost::get(&VertexPropertyTrue::phi, G_true));
      dp_t.property("z", boost::get(&VertexPropertyTrue::z, G_true));
      dp_t.property("hit_id", boost::get(&VertexPropertyTrue::hit_id, G_true));
      dp_t.property("pt_particle", boost::get(&VertexPropertyTrue::pt_particle, G_true));
      dp_t.property("eta", boost::get(&VertexPropertyTrue::eta, G_true));
      //dp_t.property("PID", boost::get(&VertexPropertyTrue::PID, G_true));
      dp_t.property("index", boost::get(&VertexPropertyTrue::index, G_true));


      std::string filenameg_t="/event"+evtid+"_TARGET.txt";
      auto path_graph_true = joinPaths(m_cfg.outputDir,  filenameg_t);
      std::fstream outGraph_TRUE(path_graph_true,std::ios::out); 
      boost::write_graphml( outGraph_TRUE, G_true,dp_t, true);
      outGraph_TRUE.close();
    }
  }    
  
  
  ACTS_INFO("end");
  return ProcessCode::SUCCESS;
}


ProcessCode GraphCreatorWriter::endRun() {

  return ProcessCode::SUCCESS;
}
