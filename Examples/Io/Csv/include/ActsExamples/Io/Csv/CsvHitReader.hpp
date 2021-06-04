// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IReader.hpp"

#include <memory>
#include <string>
#include <unordered_map>

namespace Acts {
class Surface;
}

namespace ActsExamples {

class CsvHitReader final : public IReader {
 public:
  struct Config {
    /// Where to read input files from.
    std::string inputDir;
    /// Output hits collection.
    std::string outputHits_map;
    std::string outputHits_vec;
  };

    struct Hitinformation{
        float x;
        float y;
        float z;
        uint64_t hit_id;
        uint64_t particle_id;
        Acts::GeometryIdentifier geometry_id;
        float pt_particle;
        int nhit_particle;
        int index;
        bool duplicate = false;
        bool linkUp=false;
        int nhits_wihoutDuplicate;
    };
    
  /// Construct the hit reader.
  ///
  /// @params cfg is the configuration object
  /// @params lvl is the logging level
  CsvHitReader(const Config& cfg, Acts::Logging::Level lvl);

  std::string name() const final override;

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const final override;

  /// Read out data from the input stream.
  ProcessCode read(const ActsExamples::AlgorithmContext& ctx) final override;

 private:
  Config m_cfg;
  std::pair<size_t, size_t> m_eventsRange;
  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples

