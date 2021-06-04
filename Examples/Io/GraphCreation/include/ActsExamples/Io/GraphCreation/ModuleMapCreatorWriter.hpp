// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Utilities/Logger.hpp>
#include "ActsExamples/Framework/IWriter.hpp"
#include <limits>
#include <boost/bimap.hpp>
#include <TFile.h>
#include <TTree.h>

namespace Acts {
  class TrackingVolume;
}

namespace ActsExamples {

class ModuleMapCreatorWriter : public IWriter {
 public:
  
  struct Config {
    /// The tracking geometry that is used to build the list of modules.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Where to place output files.
    std::string outputDir;
    std::string inputDir;
    std::string rootName;
    /// Cuts on the generated particles
    float minPt;
    long unsigned int minNHits;
    //flag to give cuts on geometric parameters
    bool giveCutsValues;
    
  };
  
  /// Construct the module links writer.
  ///
  /// @param cfg is the configuration object
  /// @param lvl is the logging level
  ModuleMapCreatorWriter(const Config& cfg, Acts::Logging::Level lvl);
  std::string name() const final override;
  ProcessCode write(const AlgorithmContext& context) final override;
  ProcessCode endRun() final override;

  private:
    Config m_cfg;
          
    const Acts::TrackingVolume* m_world;
    std::unique_ptr<const Acts::Logger> m_logger;

    //ModuleMap map
    typedef boost::bimap<Acts::GeometryIdentifier::Value,size_t> bimap;
    bimap m_ModuleMap;
    //Matrix of (n_module*n_module) to count numbers of possible connections
    unsigned int** m_LinkCounters;
    
    //for output
    Acts::GeometryIdentifier::Value m_tree_module1, m_tree_module2;
    unsigned int m_count;
    
    TFile *m_RootFile;
    TTree *m_treeModuleMap = new TTree("TreeModuleMap","Tree containing the module map");
    
    //adjust cuts on each variable
    float** m_module_z0_min;
    float m_z0_min;
    float** m_module_phislope_min;
    float m_phislope_min;
    float** m_module_deta_min;
    float m_deta_min;
    float** m_module_dphi_min;
    float m_dphi_min;

    float** m_module_z0_max;
    float m_z0_max;
    float** m_module_phislope_max;
    float m_phislope_max;
    float** m_module_deta_max;
    float m_deta_max;
    float** m_module_dphi_max;
    float m_dphi_max;


    void examineLayer(const Acts::TrackingVolume* vol);
          
    const Acts::Logger& logger() const { return *m_logger; }
};
}
