// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Utilities/Logger.hpp>
#include "ActsExamples/Framework/IWriter.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <limits>
#include <TFile.h>
#include <TTree.h>
#include <boost/graph/adjacency_list.hpp> // we have to test with #include <boost/graph/adjacency_matrix.hpp>
//#include <boost/property_map/dynamic_property_map.hpp>

namespace ActsExamples {

  /// Create graphs for GNN-based tracking
  ///
  class GraphCreatorWriter : public IWriter {
  public:
    struct Config {
      std::string inputDir;
      // Where to place output files.
      std::string outputDir;
      //imput of the root Module map
      std::string inputModuleMap;
      //true graph
      bool trueGraph;
      //save graph
      bool saveGraph;
      //for true flag
      float minPt;
      long unsigned int minNHits;
    };

    /// Construct the graph creator.
    /// @param cfg is the configuration object
    /// @param lvl is the logging level
    GraphCreatorWriter(const Config& cfg, Acts::Logging::Level lvl);
    std::string name() const final override;
    ProcessCode write(const AlgorithmContext& context) final override;
    ProcessCode endRun() final override;

  private:
    Config m_cfg;
    std::unique_ptr<const Acts::Logger> m_logger;                                                                                                                 
    const Acts::Logger& logger() const { return *m_logger; }
    bool compare_float(float x, float y);
    //for graphML format
    //Creation of struct to attached data to nodes and edges

    
    struct VertexProperty {
      float r;
      float phi;
      float z;
      int hit_id;
    };
    struct EdgeProperty { 
      float dEta;
      float dPhi;
      float dr;
      float dz;
    };

    //instanciate Graph
    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, VertexProperty, EdgeProperty>;
    typedef boost::graph_traits<Graph>::vertex_descriptor vertex;
    typedef boost::graph_traits<Graph>::edge_descriptor edge;

    typedef boost::unordered_map<uint64_t,int> mapB;
    
    struct VertexPropertyTrue {
      float r;
      float phi;
      float z;
      int hit_id;
      float pt_particle;
      float eta;
      int long PID;
      int index;
    };
    struct EdgePropertyTrue { 
      int is_segment_true;
      float pt_particle;
    };
    using Graph_true = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, VertexPropertyTrue,EdgePropertyTrue>;
    typedef boost::graph_traits<Graph_true>::edge_descriptor edge_true; 
    
  };

}

