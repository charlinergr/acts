// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <vector>

#include "ACTFW/Detector/IBaseDetector.hpp"
#include "ACTFW/Utilities/OptionsFwd.hpp"

namespace Acts {
class TGeoDetectorElement;
}

struct TGeoDetector : public FW::IBaseDetector {
  using DetectorElementPtr = std::shared_ptr<const Acts::TGeoDetectorElement>;
  using DetectorStore = std::vector<DetectorElementPtr>;

  /// The Store of the detector elements (lifetime: job)
  DetectorStore detectorStore;

  void addOptions(
      boost::program_options::options_description& opt) const override;

  std::pair<FW::IBaseDetector::TrackingGeometryPtr, ContextDecorators> finalize(
      const boost::program_options::variables_map& vm,
      std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) override;
};